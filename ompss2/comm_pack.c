#define ENABLE_DEBUG 0
#include "types.h"
#include "log.h"
#include "packbuf.h"
#include "dom.h"
#include "neigh.h"

#pragma oss task label("box_ghost_pack_r") \
    in(box->r) \
    in({box->neigh[i].send_rt.buf, i=0;NNEIGH}) \
    inout({box->neigh[i].send_r.buf, i=0;NNEIGH})
static void
box_ghost_pack_r(Sim *sim, Box *box)
{
    dbg("packing internal ghosts from box %d\n", box->i);

    /* Reset all PackBuf from neighbors */
    for (int i = 0; i < NNEIGH; i++) {
        packbuf_clear(&box->neigh[i].send_r);
        packbuf_debug_switch(&box->neigh[i].send_r, PB_READY, PB_PACKING);
    }

    /* Use the selection in send_rt populated by borders to pack the
     * atom position */
    for (int i = 0; i < NNEIGH; i++) {
        Neigh *neigh = &box->neigh[i];

        packbuf_debug_switch(&neigh->send_rt, PB_READY, PB_READING);
        PackBuf *pb = &neigh->send_rt;

        for (int j = 0; j < pb->natoms; j++) {
            int iatom = pb->sel[j];
            if (iatom < 0 || iatom >= box->nlocal) {
                die("atom %d outside local range\n", iatom);
            }

            Vec r = { box->r[iatom][X], box->r[iatom][Y], box->r[iatom][Z] };

            dbg("box.%d neigh.%d atom.%d (%d): sel (%e %e %e)\n",
                    box->i, neigh->i, iatom, j, r[X], r[Y], r[Z]);

            /* Enforce PBC before packing the atom position */
            if (neigh->wraps) {
                for (int d = X; d <= Z; d++)
                    r[d] += neigh->addpbc[d];

                dbg("box.%d neigh.%d atom.%d (%d): wrapped (%e %e %e)\n",
                        box->i, neigh->i, iatom, j, r[X], r[Y], r[Z]);
            }

            /* Ensure the atom is inside the destination halo domain. The atom
             * shall not exit the box domain generally, and should never exit
             * the halo domain (but it may with large time steps with no
             * re-neighboring). Only testing local boxes. TODO: test all. */
            if (ENABLE_DOMAIN_CHECK) {
                if (neigh->box) {
                    if (!in_domain(r, neigh->box->domhalo)) {
                        err("WARN: atom %d at %e %e %e"
                                " out of destination halo domain\n",
                            iatom, r[X], r[Y], r[Z]);
                    } else {
                        dbg("domain check ok\n");
                    }
                } else {
                    dbg("no domain check\n");
                }
            }

            packbuf_add(&neigh->send_r, &r, NULL, NULL);
        }

        dbg("box.%d neigh.%d: packed %d internal ghosts\n",
                box->i, neigh->i, neigh->send_r.natoms);

        packbuf_debug_switch(&neigh->send_rt, PB_READING, PB_READY);
        packbuf_debug_switch(&neigh->send_r, PB_PACKING, PB_READY);
    }
}

#pragma oss task label("neigh_border_pack_rt") \
    in(box->r) \
    out({box->neigh[i].send_rt.buf,     i=0;NNEIGH}) \
    out({box->neigh[i].send_rt.natoms,  i=0;NNEIGH}) \
    out({box->neigh[i].send_rt.sel,     i=0;NNEIGH})
static void
box_border_pack_rt(Sim *sim, Box *box)
{
    dbg("packing borders for box %2d with nlocal %d\n",
            box->i, box->nlocal);

//    for (int i = 0; i < box->nlocal; i++) {
//        if (!in_domain(box->r[i], box->domhalo)) {
//            dbg("box %d: atom %d at %e %e %e is outside halo domain\n",
//                    box->i, i, box->r[i][X], box->r[i][Y], box->r[i][Z]);
//            abort();
//        }
//    }

    /* Reset all PackBuf from neighbors */
    for (int i = 0; i < NNEIGH; i++) {
        packbuf_debug_switch(&box->neigh[i].send_rt, PB_READY, PB_PACKING);
        packbuf_clear(&box->neigh[i].send_rt);
    }

    /* Reset ghosts in this box */
    box->nghost = 0;

    for (int i = 0; i < box->nlocal; i++) {
        int delta[NDIM];

        if (in_domain_delta(box->r[i], box->domcore, delta))
            continue;

        if (ENABLE_DOMAIN_CHECK && !in_domain(box->r[i], box->dombox)) {
            die("box %d contains local atom %d at %e %e %e outside box domain\n",
                    box->i, i, box->r[i][X], box->r[i][Y], box->r[i][Z]);
        }

        Subdomain *sub = &box->sub[delta2subdom(delta)];

        for (int j = 0; j < sub->nneigh; j++) {
            Vec r = { box->r[i][X], box->r[i][Y], box->r[i][Z] };

            Neigh *neigh = sub->neigh[j];

            dbg("atom %d out of core, delta sub (%2d %2d %2d), neigh %d/%d\n",
                    i, delta[X], delta[Y], delta[Z], neigh->i,
                    sub->nneigh);

            /* Enforce PBC before packing the atom position */
            if (neigh->wraps) {
                dbg("wrapping neigh %d atom %d position from %e %e %e\n",
                        neigh->i, i, r[X], r[Y], r[Z]);
    
                for (int d = X; d <= Z; d++)
                    r[d] += neigh->addpbc[d];
    
                dbg("wrapped  neigh %d atom %d position to   %e %e %e\n",
                        neigh->i, i, r[X], r[Y], r[Z]);
            }

            /* Encode the origin of the atom in the type */
            int type = box->atomtype[i];
            packbuf_add_sel(&neigh->send_rt, &r, NULL, &type, i);
        }

    }

    for (int i = 0; i < NNEIGH; i++) {
        packbuf_debug_switch(&box->neigh[i].send_rt, PB_PACKING, PB_READY);
    }

//    if (box->i == 1) {
//        abort();
//    }

//    for (int i = 0; i < NNEIGH; i++) {
//        Neigh *neigh = &box->neigh[i];
//        dbg("box %2d neigh %2d at delta %2d %2d %2d has %8d ghosts\n",
//                box->i, neigh->i,
//                neigh->delta[X], neigh->delta[Y], neigh->delta[Z],
//                neigh->send_rt.natoms);
//    }

    dbg("rank%d: packed borders for box %2d\n", sim->rank, box->i);
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

/* Removes the atoms that lay outside the box domain and packs them in
 * the appropriate neighbor PackBuf. Only position (r), velocity (v) and
 * type (t) is copied, as the other information is not needed. Holes are
 * filled with local atoms from the end. Notice that the ghosts are
 * invalidated.*/
#pragma oss task label("box_tidy_pack_rvt") \
    inout(*(char **)&box->r) \
    out({*(char **)&box->neigh[i].send_rvt.buf, i=0;NNEIGH}) \
    out({*(char **)&box->neigh[i].send_rvt.natoms, i=0;NNEIGH})
static void
box_tidy_pack_rvt(Sim *sim, Box *box)
{
    dbg("packing out atoms for box %2d with nlocal %d\n",
            box->i, box->nlocal);

    /* Reset all PackBuf from neighbors */
    for (int i = 0; i < NNEIGH; i++) {
        packbuf_clear(&box->neigh[i].send_rvt);
        packbuf_debug_switch(&box->neigh[i].send_rvt, PB_READY, PB_PACKING);
    }

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

//        dbg("box %d packing atom %3d in neigh %d at %e %e %e\n",
//                box->i, i, neigh->i, r[X], r[Y], r[Z]);

        /* Enforce PBC before packing the atom position */
        if (neigh->wraps) {
            for (int d = X; d <= Z; d++)
                r[d] += neigh->addpbc[d];

//            dbg("               atom %3d wraps, now at %e %e %e\n",
//                    i, r[X], r[Y], r[Z]);
        }

        int type = box->atomtype[i];
        packbuf_add(&neigh->send_rvt, &r, &box->v[i], &type);

        /* Fill the hole with one atom from the end */
        int src = box->nlocal - 1, dst = i;
        //dbg("box %d moving atom %d to %d\n",
        //        box->i, src, dst);
        copy_atom_rvt(box, src, dst);
        box->nlocal--;

    }

    for (int i = 0; i < NNEIGH; i++) {
        packbuf_debug_switch(&box->neigh[i].send_rvt, PB_PACKING, PB_READY);
    }

    for (int i = 0; i < NNEIGH; i++) {
        Neigh *neigh = &box->neigh[i];
        if (neigh->send_rvt.natoms > 0) {
            dbg("packed %d atoms in box%d:neigh%02d pb=%p\n",
                    neigh->send_rvt.natoms,
                    box->i, neigh->i, &neigh->send_rvt);
        }
    }

    dbg("after packing out atoms for box %2d, nlocal=%d\n",
            box->i, box->nlocal);
}

void
comm_pack(Sim *sim, enum pb_type type)
{
    dbg("comm_pack %s\n", PB_TYPENAME(type));

    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];

        switch (type) {
        case PB_SEND_R: box_ghost_pack_r(sim, box); break;
        case PB_SEND_RT: box_border_pack_rt(sim, box); break;
        case PB_SEND_RVT: box_tidy_pack_rvt(sim, box); break;
        default: die("not implemented\n");
        }
    }

    #pragma oss taskwait
}
