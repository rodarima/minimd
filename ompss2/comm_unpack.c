#define ENABLE_DEBUG 1
#include "types.h"
#include "log.h"
#include "packbuf.h"
#include "box.h"
#include "dom.h"
#include "neigh.h"

#include <math.h>

static void
neigh_ghost_unpack_r(Sim *sim, Box *box, Neigh *neigh)
{
    packbuf_debug_switch(&neigh->recv_r, PB_READY, PB_RECVING);
    //packbuf_debug_switch(&neigh->recv_rt, PB_READY, PB_READING);

    if (neigh->recv_rt.natoms != neigh->recv_r.natoms)
        abort();

    if (neigh->recv_r.natoms != 0) {
        int nnew = neigh->recv_r.natoms;
        int ntot = box->nlocal + box->nghost + nnew;
        if (ntot > box->nalloc) {
            die("error: box%d:neigh%d cannot unpack %d atoms, capacity exeeeded\n",
                    box->i, neigh->i, nnew);
        }

        int i = box->nlocal + box->nghost;

        /* The unpack order must be kept the same to match the ghost atom
         * order given by borders */
        packbuf_unpack(&neigh->recv_r, &box->r[i], NULL, NULL);

        /* We cannot check the domain bounds of the new ghost atom
         * positions, as they are moving around, even exceeding the halo
         * domain */

        box->nghost += neigh->recv_r.natoms;
    }

    packbuf_debug_switch(&neigh->recv_r, PB_RECVING, PB_READY);
    //packbuf_debug_switch(&neigh->recv_rt, PB_READING, PB_READY);
}

#pragma oss task label("box_ghost_unpack_r_neigh") \
    in({box->neigh[i].recv_rt.natoms,    i=0;NNEIGH}) \
    in({box->neigh[i].recv_rt.buf,       i=0;NNEIGH}) \
    in({box->neigh[i].recv_r.natoms,     i=0;NNEIGH}) \
    in({box->neigh[i].recv_r.buf,        i=0;NNEIGH}) \
    inout(box->r)
static void
box_ghost_unpack_r(Sim *sim, Box *box)
{
    dbg("unpacking r for box %d\n", box->i);

    int old_nghost = box->nghost;
    box->nghost = 0;

    for (int j = 0; j < NNEIGH; j++) {
        /* Use the inverse order for unpack */
        Neigh *opp = box->neigh[j].opposite;
        neigh_ghost_unpack_r(sim, box, opp);
    }

    if (ENABLE_ATOM_COUNT_CHECK) {
        /* Wait until all unpack have finished */
//        #pragma oss task label("box_ghost_unpack_r:atomcheck") \
//            in(*(char **)&box->r)
        if (box->nghost != old_nghost) {
            die("nghost atoms don't match %d != %d\n",
                    box->nghost, old_nghost);
        }
    }
}

/* FIXME: We shouldn't need to use inout */
#pragma oss task label("neigh_border_unpack_rt") \
    inout(neigh->recv_rt.buf) \
    inout(neigh->recv_rt.natoms) \
    out(neigh->recv_r.recvnatoms) \
    out(box->r)
static void
neigh_border_unpack_rt(Sim *sim, Box *box, Neigh *neigh)
{
    neigh->recv_r.recvnatoms = neigh->recv_rt.natoms;

    if (neigh->recv_rt.natoms == 0)
        return;

    packbuf_debug_switch(&neigh->recv_rt, PB_READY, PB_UNPACKING);

    /* Ensure we have room to place the new ghost atoms */
    int nnew = neigh->recv_rt.natoms;
    int nend = box->nlocal + box->nghost;
    int ntot = nend + nnew;
    int oldalloc = box->nalloc;
    box_realloc(box, ntot);

//    dbg("unpacking %d atoms from neigh %d into box %d (%d -> %d)\n",
//            neigh->recv_rt.natoms, neigh->i, box->i, nend, ntot);

    //dbg("box %d realloc from %d to %d (nnew=%d nend=%d ntot=%d)\n",
    //        box->i, oldalloc, box->nalloc, nnew, nend, ntot);

    /* Unpack the position and type at the end of the local atoms */
    packbuf_unpack(&neigh->recv_rt, &box->r[nend], NULL, &box->atomtype[nend]);

    for (int i = nend; i < ntot; i++) {
        Vec r = { box->r[i][X], box->r[i][Y], box->r[i][Z] };
        if (ENABLE_DOMAIN_CHECK) {
            if (!in_domain(r, box->domhalo)) {
                die("rank%d.box%d.neigh%d: unpacked ghost atom %d at %e %e %e outside halo domain\n",
                    sim->rank, box->i, neigh->i, i, r[X], r[Y], r[Z]);
            }
            if (in_domain(r, box->dombox)) {
                die("rank%d.box%d.neigh%d: unpacked ghost atom %d at %e %e %e inside box domain\n",
                        sim->rank, box->i, neigh->i, i, r[X], r[Y], r[Z]);
            }
        }
    }

    /* Adjust the number of ghost atoms in the box */
    box->nghost += nnew;

    packbuf_debug_switch(&neigh->recv_rt, PB_UNPACKING, PB_READY);
}

static void
box_border_unpack_rt(Sim *sim, Box *box)
{
    dbg("unpacking rt for box %d\n", box->i);

    for (int j = 0; j < NNEIGH; j++) {
        Neigh *neigh = &box->neigh[j];
        neigh_border_unpack_rt(sim, box, neigh);
    }
}

static void
check_atom(Sim *sim, Box *box, Vec r)
{
    dbg("matching incoming atom at %e %e %e\n",
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

        if (jindbin < 0 || jindbin >= box->nbinsalloc)
            die("near atom bin is outside the range\n");

        Bin *jbin = &box->bin[jindbin];

        /* Find atoms within 0.3 R_neigh in the bin */
        for (int k = 0; k < jbin->natoms; k++) {
            int jatom = jbin->atom[k];

            Vec rj = { box->r[jatom][X], box->r[jatom][Y], box->r[jatom][Z] };

            double dist_sq = get_distsq(r, rj);
            double dist = sqrt(dist_sq);

            if (dist < closest)
                closest = dist;

            if (dist < closestghost && jatom >= box->nlocal)
                closestghost = dist;

            /* Ignore atoms which are not too close */
            if (dist >= limit) {
                dbg("ignoring atom %d in bin %d with dist %f\n",
                        jatom, jindbin, dist);
                continue;
            }

            if (jatom < box->nlocal) {
                err("WARNING: too close to local %d at %e %e %e\n",
                        jatom, rj[X], rj[Y], rj[Z]);
            } else {
                dbg("potential match with ghost %d with dist %e\n",
                        jatom, dist);
                match++;
            }
        }
    }

    if (match != 1) {
        die("no match, closest %f (ghost %f), limit %f\n",
                closest, closestghost, limit);
    }
}


#pragma oss task label("neigh_tidy_unpack_rvt") \
    inout(*(char **)&box->r) \
    inout(*(char **)&box->v) \
    inout(*(char **)&box->f) /* May realloc f too */\
    inout(*(char **)&neigh->recv_rvt.buf) \
    inout(*(char **)&neigh->recv_rvt.natoms)
static void
neigh_tidy_unpack_rvt(Sim *sim, Box *box, Neigh *neigh)
{
    packbuf_debug_switch(&neigh->recv_rvt, PB_READY, PB_UNPACKING);
//    if (neigh->recv_rvt.natoms > 0) {
//        dbg("box %d: unpacking %d atoms from neigh %d\n",
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
    packbuf_debug_switch(&neigh->recv_rvt, PB_UNPACKING, PB_READY);
}

static void
box_tidy_unpack_rvt(Sim *sim, Box *box)
{
    dbg("unpacking rvt for box %d\n", box->i);

    for (int j = 0; j < NNEIGH; j++) {
        Neigh *neigh = &box->neigh[j];
        neigh_tidy_unpack_rvt(sim, box, neigh);
    }
}

void
comm_unpack(Sim *sim, enum pb_type type)
{
    dbg("comm_unpack %s\n", PB_TYPENAME(type));

    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];

        switch (type) {
        case PB_RECV_R: box_ghost_unpack_r(sim, box); break;
        case PB_RECV_RT: box_border_unpack_rt(sim, box); break;
        case PB_RECV_RVT: box_tidy_unpack_rvt(sim, box); break;
        default: die("not implemented\n");
        }
    }

    #pragma oss taskwait
}
