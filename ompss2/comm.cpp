#include "types.h"
#include "neigh.h"
#include "dom.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

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
//    fprintf(stderr, "packing out atoms for box %2d with nlocal %d\n",
//            box->i, box->nlocal);

    /* Reset all PackBuf from neighbors */
    for (int i = 0; i < NNEIGH; i++)
        packbuf_clear(&box->neigh[i].send_rvt);

    /* Invalidate ghosts, as we are going to modify box->nlocal */
    box->nghost = -666;

    for (int i = 0; i < box->nlocal; ) {

        Vec r = { box->r[i][X], box->r[i][Y], box->r[i][Z] };
        int delta[NDIM];

        if (in_domain_delta(r, box->dombox, delta)) {
            i++;
            continue;
        }

        Neigh *neigh = &box->neigh[delta2neigh(delta)];

//        fprintf(stderr,     "box %d packing atom %3d in neigh %d at %e %e %e\n",
//                box->i, i, neigh->i, r[X], r[Y], r[Z]);

        /* Enforce PBC before packing the atom position */
        if (neigh->wraps) {
            for (int d = X; d <= Z; d++)
                r[d] += neigh->addpbc[d];

//            fprintf(stderr, "               atom %3d wraps, now at %e %e %e\n",
//                    i, r[X], r[Y], r[Z]);
        }

        int type = box->atomtype[i] + 10000 * (box->i + 50);
        packbuf_add(&neigh->send_rvt, &r, &box->v[i], &type);

        /* Fill the hole with one atom from the end */
        int src = box->nlocal - 1, dst = i;
        //fprintf(stderr, "box %d moving atom %d to %d\n",
        //        box->i, src, dst);
        copy_atom_rvt(box, src, dst);
        box->nlocal--;

    }

//    for (int i = 0; i < NNEIGH; i++) {
//        Neigh *neigh = &box->neigh[i];
//        if (neigh->send_rvt.natoms > 0) {
//            fprintf(stderr, "packed %d atoms in box%d:neigh%02d\n",
//                    neigh->send_rvt.natoms,
//                    box->i, neigh->i);
//        }
//    }
}

static void
box_send_rvt(Sim *sim, Box *box, Neigh *neigh)
{
//    if (neigh->send_rvt.natoms > 0) {
//        fprintf(stderr, "box %d sending %d atoms to neigh %d\n",
//                box->i, neigh->send_rvt.natoms, neigh->i);
//    }

    if (neigh->rank != sim->rank) {
        /* Use MPI for inter process comm */
        packbuf_mpisend(&neigh->send_rvt, neigh->rank, neigh->i);
    } else {
        /* No-op: will be copied in box_recv_rvt */
    }
}

static void
box_recv_rvt(Sim *sim, Box *dstbox, Neigh *dstneigh)
{
    /* Receive */
    if (dstneigh->rank != sim->rank) {
        /* Use MPI for inter process comm */
        packbuf_mpirecv(&dstneigh->recv_rvt, dstneigh->rank, dstneigh->i);
    } else {
        /* Shared memory for intra-process. This can be avoided if
         * we pack directly into the receiving buffer. */
        Box *srcbox = dstneigh->box;
        Neigh *srcneigh = &srcbox->neigh[opposite_neigh(dstneigh->i)];
        packbuf_shmcopy(&srcneigh->send_rvt, &dstneigh->recv_rvt);

//        if (srcneigh->send_rvt.natoms > 0) {
//            fprintf(stderr, "shmcopy %d atoms from box%d:neigh%02d -> box%d:neigh%02d\n",
//                    srcneigh->send_rvt.natoms,
//                    srcbox->i, srcneigh->i,
//                    dstbox->i, dstneigh->i);
//        }
    }

}

static void
check_atom(Sim *sim, Box *box, Vec r)
{
    fprintf(stderr, "matching incoming atom at %e %e %e\n",
            r[X], r[Y], r[Z]);

    if(!in_domain(r, box->dombox))
        abort();

    /* Identify the atom bin */
    int iindbin = get_atom_bin(sim, box, r);
    Bin *ibin = &box->bin[iindbin];


    int match = 0;
    double closest = 1e50;
    double closestghost = 1e50;
    double limit = 0.1 * sim->R_neigh;
    /* Find nearing atoms in the stencil bins */
    for (int j = 0; j < box->nstencil; j++) {
        int jindbin = iindbin + box->stencil[j];

        if (jindbin < 0 || jindbin >= box->nbinsalloc) {
            fprintf(stderr, "near atom bin is outside the range\n");
            abort();
        }

        Bin *jbin = &box->bin[jindbin];

        /* Find atoms within 0.3 R_neigh in the bin */
        for (int k = 0; k < jbin->natoms; k++) {
            int jatom = jbin->atom[k];

            Vec rj = { box->r[jatom][X], box->r[jatom][Y], box->r[jatom][Z] };
            int jtype = box->atomtype[jatom];

            double dist_sq = get_distsq(r, rj);
            double dist = sqrt(dist_sq);

            if (dist < closest)
                closest = dist;

            if (dist < closestghost && jatom >= box->nlocal)
                closestghost = dist;

            /* Ignore atoms which are not too close */
            if (dist >= limit) {
                fprintf(stderr, "ignoring atom %d in bin %d with dist %f\n",
                        jatom, jindbin, dist);
                continue;
            }

            if (jatom < box->nlocal) {
                fprintf(stderr, "WARNING: too close to local %d at %e %e %e\n",
                        jatom, rj[X], rj[Y], rj[Z]);
            } else {
                fprintf(stderr, "potential match with ghost %d with dist %e\n",
                        jatom, dist);
                match++;
            }
        }
    }

    if (match != 1) {
        fprintf(stderr, "no match, closest %f (ghost %f), limit %f\n",
                closest, closestghost, limit);
        abort();
    }
}

static void
box_unpack_rvt(Sim *sim, Box *box, Neigh *neigh)
{
//    if (neigh->recv_rvt.natoms > 0) {
//        fprintf(stderr, "box %d: unpacking %d atoms from neigh %d\n",
//                box->i, neigh->recv_rvt.natoms, neigh->i);
//    }
    /* Ensure we have room to place the new local atoms */
    int n = box->nlocal + neigh->recv_rvt.natoms;
    box_realloc(box, n);

//    /* Here we are receiving new local atoms that have just moved into
//     * our box. Ensure that there was a ghost atom before and that is
//     * not too close to an already existing atom */
//    for (int i = 0; i < neigh->recv_rvt.natoms; i++) {
//        double *r = &neigh->recv_rvt.buf[i * (NDIM * 2 + 1)];
//        check_atom(sim, box, *(Vec *) r);
//    }

    /* Unpack at the end of the local atoms */
    Vec *r = &box->r[box->nlocal];
    Vec *v = &box->v[box->nlocal];
    int *types = &box->atomtype[box->nlocal];
    packbuf_unpack(&neigh->recv_rvt, r, v, types);

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
        }
    }

    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int i = 0; i < NNEIGH; i++) {
            Neigh *neigh = &box->neigh[i];
            box_recv_rvt(sim, box, neigh->opposite);
        }
    }

    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int i = 0; i < NNEIGH; i++) {
            Neigh *neigh = &box->neigh[i];
            box_unpack_rvt(sim, box, neigh->opposite);
        }
    }
}

static void
box_pack_borders(Sim *sim, Box *box)
{
//    fprintf(stderr, "packing borders for box %2d with nlocal %d\n",
//            box->i, box->nlocal);

//    for (int i = 0; i < box->nlocal; i++) {
//        if (!in_domain(box->r[i], box->domhalo)) {
//            fprintf(stderr, "box %d: atom %d at %e %e %e is outside halo domain\n",
//                    box->i, i, box->r[i][X], box->r[i][Y], box->r[i][Z]);
//            abort();
//        }
//    }

    /* Reset all PackBuf from neighbors */
    for (int i = 0; i < NNEIGH; i++) {
        packbuf_clear(&box->neigh[i].send_rt);
    }

    /* Reset ghosts in this box */
    box->nghost = 0;

    for (int i = 0; i < box->nlocal; i++) {
        int delta[NDIM];

        if (in_domain_delta(box->r[i], box->domcore, delta))
            continue;

        if (ENABLE_DOMAIN_CHECK && !in_domain(box->r[i], box->dombox)) {
            fprintf(stderr, "box %d contains local atom %d at %e %e %e outside box domain\n",
                    box->i, i, box->r[i][X], box->r[i][Y], box->r[i][Z]);
            abort();
        }

        Subdomain *sub = &box->sub[delta2subdom(delta)];

        for (int j = 0; j < sub->nneigh; j++) {
            Vec r = { box->r[i][X], box->r[i][Y], box->r[i][Z] };

            Neigh *neigh = sub->neigh[j];

//            fprintf(stderr, "atom %d out of core, delta sub (%2d %2d %2d), neigh %d/%d\n",
//                    i, delta[X], delta[Y], delta[Z], neigh->i,
//                    sub->nneigh);

            /* Enforce PBC before packing the atom position */
            if (neigh->wraps) {
//                fprintf(stderr, "wrapping neigh %d atom %d position from %e %e %e\n",
//                        neigh->i, i, r[X], r[Y], r[Z]);
    
                for (int d = X; d <= Z; d++)
                    r[d] += neigh->addpbc[d];
    
//                fprintf(stderr, "wrapped  neigh %d atom %d position to   %e %e %e\n",
//                        neigh->i, i, r[X], r[Y], r[Z]);
            }

            /* Encode the origin of the atom in the type */
            int type = box->atomtype[i] + (box->i + 1) * 10000;
            packbuf_add_sel(&neigh->send_rt, &r, NULL, &type, i);
        }

    }
//    if (box->i == 1) {
//        abort();
//    }

//    for (int i = 0; i < NNEIGH; i++) {
//        Neigh *neigh = &box->neigh[i];
//        fprintf(stderr, "box %2d neigh %2d at delta %2d %2d %2d has %8d ghosts\n",
//                box->i, neigh->i,
//                neigh->delta[X], neigh->delta[Y], neigh->delta[Z],
//                neigh->send_rt.natoms);
//    }
}

static void
box_unpack_rt(Sim *sim, Box *box, Neigh *neigh)
{
    if (neigh->recv_rt.natoms == 0)
        return;

    /* Ensure we have room to place the new ghost atoms */
    int nnew = neigh->recv_rt.natoms;
    int nend = box->nlocal + box->nghost;
    int ntot = nend + nnew;
    int oldalloc = box->nalloc;
    box_realloc(box, ntot);

//    fprintf(stderr, "unpacking %d atoms from neigh %d into box %d (%d -> %d)\n",
//            neigh->recv_rt.natoms, neigh->i, box->i, nend, ntot);

    //fprintf(stderr, "box %d realloc from %d to %d (nnew=%d nend=%d ntot=%d)\n",
    //        box->i, oldalloc, box->nalloc, nnew, nend, ntot);

    /* Unpack the position and type at the end of the local atoms */
    packbuf_unpack(&neigh->recv_rt, &box->r[nend], NULL, &box->atomtype[nend]);

    for (int i = nend; i < ntot; i++) {
        Vec r = { box->r[i][X], box->r[i][Y], box->r[i][Z] };
        /* Relaxed */
//        if (!in_domain(r, box->domhalo)) {
//            fprintf(stderr, "error: unpacked ghost atom %d at %e %e %e outside halo domain\n",
//                    i, r[X], r[Y], r[Z]);
//            abort();
//        }
        if (ENABLE_DOMAIN_CHECK && in_domain(r, box->dombox)) {
            fprintf(stderr, "error: unpacked ghost atom %d at %e %e %e inside box domain\n",
                    i, r[X], r[Y], r[Z]);
            abort();
        }
    }

    /* Adjust the number of ghost atoms in the box */
    box->nghost += nnew;
}

static void
comm_borders_pack(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];

        box_pack_borders(sim, box);
    }
}

static void
comm_borders_send(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            Neigh *neigh = &box->neigh[j];
            if (neigh->rank != sim->rank) {
                packbuf_mpisend(&neigh->send_rt, neigh->rank, neigh->i);
            } else {
                /* No-op: will be copied via shared memory at recv */
            }
        }
    }
}

static void
comm_borders_recv(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            /* Use the inverse order for reception */
            Neigh *dstneigh = box->neigh[j].opposite;

            /* Clear receive buffer */
            packbuf_clear(&dstneigh->recv_rt);

            if (dstneigh->rank != sim->rank) {
                packbuf_mpirecv(&dstneigh->recv_rt, dstneigh->rank, dstneigh->i);
            } else {
                /* Get the source box from the neighbor and find the
                 * neighbor which contains the send buffer. Example:
                 *
                 * +Y
                 * ^
                 * |   +-----+         In this example, box_0 must receive
                 * |   |box_1|         the atoms from box_1. The neighbor
                 * |   |     |         neigh_15 pointing upwards (+Y) is
                 *     +--|--+         where the recv buffer is located.
                 *        v neigh_10   
                 *          (0, -1, 0) However, we must first access the
                 *                     neigh_15->box to find box_1, and
                 *          neigh_15   then compute the opposite neighbor.
                 *        ^ (0, +1, 0) 
                 *     +--|--+         The opposite neighbor is at
                 *     |box_0|         neigh_10 = opposite_neigh(15)
                 *     |     |         The send buffer is then used from
                 *     +-----+         neigh_10.
                 *
                 */
                Box *srcbox = dstneigh->box;
                Neigh *srcneigh = &srcbox->neigh[opposite_neigh(dstneigh->i)];
                packbuf_shmcopy(&srcneigh->send_rt, &dstneigh->recv_rt);
            }
        }
    }
}

static void
comm_borders_unpack(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            /* Use the inverse order for unpack */
            Neigh *opp = box->neigh[j].opposite;
            box_unpack_rt(sim, box, opp);
        }

        if (box->nghost == 0)
            abort();
    }
}

void
comm_borders(Sim *sim)
{
    comm_borders_pack(sim);
    comm_borders_send(sim);
    comm_borders_recv(sim);
    comm_borders_unpack(sim);
}

static void
box_pack_ghost_r(Sim *sim, Box *box)
{
//    fprintf(stderr, "packing internal ghosts from box %d\n", box->i);

    /* Reset all PackBuf from neighbors */
    for (int i = 0; i < NNEIGH; i++) {
        packbuf_clear(&box->neigh[i].send_r);
    }

    /* Use the selection in send_rt populated by borders to pack the
     * atom position */
    for (int i = 0; i < NNEIGH; i++) {
        Neigh *neigh = &box->neigh[i];
        PackBuf *pb = &neigh->send_rt;

        for (int j = 0; j < pb->natoms; j++) {
            int iatom = pb->sel[j];
            if (iatom < 0 || iatom >= box->nlocal) {
                fprintf(stderr, "atom %d outside local range\n", iatom);
                abort();
            }

            Vec r = { box->r[iatom][X], box->r[iatom][Y], box->r[iatom][Z] };

            /* Enforce PBC before packing the atom position */
            if (neigh->wraps) {
                for (int d = X; d <= Z; d++)
                    r[d] += neigh->addpbc[d];
            }

            packbuf_add(&neigh->send_r, &r, NULL, NULL);
        }
    }
}

static void
send_ghost_position(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            Neigh *neigh = &box->neigh[j];
            if (neigh->rank != sim->rank) {
                /* Only send the atom buffer, the receive end already
                 * knows the size */
                packbuf_mpisend_buf(&neigh->send_r, neigh->rank, neigh->i);
            } else {
                /* No-op: will be copied via shared memory at recv */
            }
        }
    }
}

static void
recv_ghost_position(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            /* Use the inverse order for reception */
            Neigh *dstneigh = box->neigh[j].opposite;

            /* Get the number of atoms to be received from the pack
             * buffer used in the borders */
            int natoms = dstneigh->recv_rt.natoms;

            /* Clear receive buffer */
            packbuf_clear(&dstneigh->recv_r);

            if (dstneigh->rank != sim->rank) {
                packbuf_mpirecv_buf(&dstneigh->recv_r, dstneigh->rank,
                        dstneigh->i, natoms);
            } else {
                /* Get the source box from the neighbor and find the
                 * neighbor which contains the send buffer */
                Box *srcbox = dstneigh->box;
                Neigh *srcneigh = &srcbox->neigh[opposite_neigh(dstneigh->i)];

                if (srcneigh->send_r.natoms != natoms)
                    abort();

                packbuf_shmcopy(&srcneigh->send_r, &dstneigh->recv_r);
            }

//            if (natoms != 0) {
//                fprintf(stderr, "box %d: unpacked %d ghost atoms from neigh %d\n",
//                        box->i, natoms, j);
//            }
        }
    }
}

static void
box_unpack_r(Sim *sim, Box *box, Neigh *neigh)
{
    if (neigh->recv_r.natoms == 0)
        return;

    if (neigh->recv_rt.natoms != neigh->recv_r.natoms)
        abort();

    int i = box->nlocal + box->nghost;

    /* The unpack order must be kept the same to match the ghost atom
     * order given by borders */
    packbuf_unpack(&neigh->recv_r, &box->r[i], NULL, NULL);

    /* We cannot check the domain bounds of the new ghost atom
     * positions, as they are moving around, even exceeding the halo
     * domain */

    box->nghost += neigh->recv_r.natoms;
}

static void
unpack_ghost_position(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];

        int old_nghost = box->nghost;
        box->nghost = 0;

        for (int j = 0; j < NNEIGH; j++) {
            /* Use the inverse order for unpack */
            Neigh *opp = box->neigh[j].opposite;
            box_unpack_r(sim, box, opp);
        }

        if (box->nghost != old_nghost) {
            fprintf(stderr, "nghost atoms don't match %d != %d\n",
                    box->nghost, old_nghost);
            abort();
        }
    }
}

/* Send/recv ghost positions (communicate) */
void
comm_ghost_position(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        box_pack_ghost_r(sim, box);
    }

    send_ghost_position(sim);
    recv_ghost_position(sim);

    unpack_ghost_position(sim);
}
