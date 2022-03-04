#include "types.h"
#include "neigh.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#ifdef USE_TAMPI
#include <TAMPI.h>
#else
int TAMPI_Iwait(MPI_Request *r, MPI_Status *s) {
  return MPI_Wait(r, s);
}

int TAMPI_Iwaitall(int n, MPI_Request *r, MPI_Status *s) {
  return MPI_Waitall(n, r, s);
}
#endif


/* Ensures that at least n atoms fit in the box (both local and ghosts) */
static void
box_realloc(Box *box, int n)
{
    if (box->nalloc >= n)
        return;

    box->r = (Vec *) safe_realloc(box->r, n * sizeof(box->r[0]));
    box->v = (Vec *) safe_realloc(box->v, n * sizeof(box->v[0]));
    box->f = (Vec *) safe_realloc(box->f, n * sizeof(box->f[0]));
    box->atomtype = (int *) safe_realloc(box->atomtype, n * sizeof(box->atomtype[0]));
    box->nalloc = n;
}

static void
copy_atom_rvt(Box *box, int src, int dst)
{
    for (int d = X; d <= Z; d++)
        box->r[dst][d] = box->r[src][d];

    for (int d = X; d <= Z; d++)
        box->v[dst][d] = box->v[src][d];

    box->atomtype[dst] = box->atomtype[src];
}

/* Returns 0 if the position is ouside the domain, and writes the offset
 * in delta. Otherwise returns 1 */
static int
is_inside_domain(Vec r, Domain dom, int delta[NDIM])
{
    int inside = 1;

    for (int d = X; d <= Z; d++) {
        if (r[d] < dom[d][LO]) {
            inside = 0;
            delta[d] = -1;
        } else if (r[d] >= dom[d][HI]) {
            inside = 0;
            delta[d] = +1;
        } else {
            delta[d] = 0;
        }
    }

    return inside;
}

static int
build_tag(Box *box, Neigh *neigh)
{
    int tag = box->i * NNEIGH + neigh->i;

    /* Ensure the tag is within the MPI standard limit */
    if (tag >= 32767) {
        fprintf(stderr, "tag exceed limit: %d >= %d\n", tag, 32767);
        abort();
    }

    return tag;
}

/* Removes the atoms that lay outside the box domain and packs them in
 * the appropriate neighbor PackBuf. Only position (r), velocity (v) and
 * type (t) is copied, as the other information is not needed. Holes are
 * filled with local atoms from the end. Notice that the ghosts are
 * invalidated.*/
static void
box_pack_rvt(Sim *sim, Box *box)
{
    fprintf(stderr, "packing out atoms for box %2d with nlocal %d\n",
            box->i, box->nlocal);

    /* Reset all PackBuf from neighbors */
    for (int i = 0; i < NNEIGH; i++)
        packbuf_clear(&box->neigh[i].send_rvt);

    /* Invalidate ghosts, as we are going to modify box->nlocal */
    box->nghost = -666;

    for (int i = 0; i < box->nlocal; ) {

        Vec r = { box->r[i][X], box->r[i][Y], box->r[i][Z] };
        int delta[NDIM];

        if (is_inside_domain(r, box->dombox, delta)) {
            i++;
            continue;
        }

        //fprintf(stderr, "box %d packing atom %d\n", box->i, i);

        Neigh *neigh = &box->neigh[delta2neigh(delta)];

        /* Enforce PBC before packing the atom position */
        if (neigh->wraps) {
            for (int d = X; d <= Z; d++)
                r[d] += neigh->addpbc[d];
        }

        packbuf_add_rvt(&neigh->send_rvt, r, box->v[i], box->atomtype[i]);

        int src = box->nlocal - 1, dst = i;
        copy_atom_rvt(box, src, dst);
        box->nlocal--;
    }
}

static void
box_send_rvt(Sim *sim, Box *box, Neigh *neigh)
{
    fprintf(stderr, "box %d sending %d atoms to neigh %d with tag %d\n",
            box->i, neigh->send_rvt.natoms, neigh->i, neigh->i);

    if (neigh->rank != sim->rank) {
        /* Use MPI for inter process comm */
        packbuf_mpisend(&neigh->send_rvt, neigh->rank, neigh->i);
    } else {
        /* Shared memory for intra-process. This can be avoided if
         * we pack directly into the receiving buffer. */
        packbuf_shmcopy(&neigh->send_rvt, &neigh->opposite->recv_rvt);
    }
}

static void
box_recv_rvt(Sim *sim, Box *box, Neigh *neigh)
{
    /* Receive */
    if (neigh->rank != sim->rank) {
        /* Use MPI for inter process comm */
        packbuf_mpirecv(&neigh->recv_rvt, neigh->rank, neigh->i);
    } else {
        /* No-op: already in recv_rvt */
    }
    fprintf(stderr, "box %d received %d atoms from neigh %d with tag %d\n",
            box->i, neigh->recv_rvt.natoms, neigh->i, neigh->i);
}

static void
box_unpack_rvt(Sim *sim, Box *box, Neigh *neigh)
{
    /* Ensure we have room to place the new local atoms */
    int n = box->nlocal + neigh->recv_rvt.natoms;
    box_realloc(box, n);

    /* Unpack at the end of the local atoms */
    Vec *r = &box->r[box->nlocal];
    Vec *v = &box->v[box->nlocal];
    int *types = &box->atomtype[box->nlocal];
    packbuf_unpack_rvt(&neigh->recv_rvt, r, v, types);

    /* Adjust the number of local atoms in the box */
    box->nlocal = n;
}

void
comm_atoms_correct_box(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        box_pack_rvt(sim, box);
    }

    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int i = 0; i < NNEIGH; i++) {
            Neigh *neigh = &box->neigh[i];

            box_send_rvt(sim, box, neigh);
            box_recv_rvt(sim, box, neigh->opposite);
            box_unpack_rvt(sim, box, neigh->opposite);
        }
    }
}

static void
box_pack_borders(Sim *sim, Box *box)
{
    fprintf(stderr, "packing borders for box %2d with nlocal %d\n",
            box->i, box->nlocal);

    /* Reset all PackBuf from neighbors */
    for (int i = 0; i < NNEIGH; i++) {
        packbuf_clear(&box->neigh[i].send_rt);
        packbuf_clear(&box->neigh[i].recv_rt);
    }

    /* Reset ghosts in this box */
    box->nghost = 0;

    for (int i = 0; i < box->nlocal; i++) {
        int delta[NDIM];

        if (is_inside_domain(box->r[i], box->domcore, delta))
            continue;

        //fprintf(stderr, "atom %d out of core, delta (%2d %2d %2d)\n",
        //        i, delta[X], delta[Y], delta[Z]);

        Neigh *neigh = &box->neigh[delta2neigh(delta)];

        //fprintf(stderr, "atom %d in neigh %d\n",
        //        i, neigh->i);

        packbuf_add_rt(&neigh->send_rt, box->r[i], box->atomtype[i]);
    }

    for (int i = 0; i < NNEIGH; i++) {
        Neigh *neigh = &box->neigh[i];
        fprintf(stderr, "box %2d neigh %2d at delta %2d %2d %2d has %8d ghosts\n",
                box->i, neigh->i,
                neigh->delta[X], neigh->delta[Y], neigh->delta[Z],
                neigh->send_rt.natoms);
    }
}

static void
box_unpack_rt(Sim *sim, Box *box, Neigh *neigh)
{
    /* Ensure we have room to place the new ghost atoms */
    int nnew = neigh->recv_rt.natoms;
    int nend = box->nlocal + box->nghost;
    int ntot = nend + nnew;
    int oldalloc = box->nalloc;
    box_realloc(box, ntot);

    fprintf(stderr, "unpacking %d atoms from neigh %d into box %d (%d -> %d)\n",
            neigh->recv_rt.natoms, neigh->i, box->i, nend, ntot);

    //fprintf(stderr, "box %d realloc from %d to %d (nnew=%d nend=%d ntot=%d)\n",
    //        box->i, oldalloc, box->nalloc, nnew, nend, ntot);

    /* Unpack the position and type at the end of the local atoms */
    Vec *r = &box->r[nend];
    int *types = &box->atomtype[nend];
    packbuf_unpack_rt(&neigh->recv_rt, r, types);

    /* Adjust the number of ghost atoms in the box */
    box->nghost += nnew;
}

void
comm_borders(Sim *sim)
{
    fprintf(stderr, "comm_borders begins\n");
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        box_pack_borders(sim, box);
    }

    for (int i = 0; i < sim->nboxes; i++) {
        fprintf(stderr, "comm_borders processing box %d\n", i);
        Box *box = &sim->box[i];
        for (int i = 0; i < NNEIGH; i++) {
            Neigh *neigh = &box->neigh[i];
            Neigh *opp = neigh->opposite;

            if (neigh->rank != sim->rank) {
                packbuf_mpisend(&neigh->send_rt, neigh->rank, neigh->i);
                packbuf_mpirecv(&opp->recv_rt, opp->rank, opp->i);
            } else {
                fprintf(stderr, "copying %d atoms from neigh %d to %d\n",
                        neigh->send_rt.natoms, neigh->i, opp->i);

                packbuf_shmcopy(&neigh->send_rt, &opp->recv_rt);
            }

            box_unpack_rt(sim, box, opp);
        }
    }
    fprintf(stderr, "comm_borders ends\n");
}
