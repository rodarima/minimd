#define ENABLE_DEBUG 0
#include "comm.h"
#include "types.h"
#include "log.h"
#include "packbuf.h"
#include "dom.h"
#include "neigh.h"
#include "box.h"

#pragma oss task label("box_ghost_pack_r") \
    in(box->r) \
    inout({box->pb[PB_RT][PB_SEND][i]->sel,     i=0;NNEIGH}) \
    inout({box->pb[PB_RT][PB_SEND][i]->natoms,  i=0;NNEIGH}) \
    inout({box->pb[PB_R ][PB_SEND][i]->data,     i=0;NNEIGH})
static void
box_ghost_pack_r(Sim *sim, Box *box)
{
    dbg("packing internal ghosts from box %d\n", box->i);

    /* Reset all PackBuf from neighbors */
    box_packbuf_clear(box, PB_R, PB_SEND);
    box_packbuf_switch(box, PB_R, PB_SEND, PB_READY, PB_PACKING);
    box_packbuf_switch(box, PB_RT, PB_SEND, PB_READY, PB_READING);

    for (int i = 0; i < NNEIGH; i++) {
        Neigh *neigh = &box->neigh[i];
        PackBuf *pb_r = box->pb[PB_R][PB_SEND][i];
        PackBuf *pb_rt = box->pb[PB_RT][PB_SEND][i];

        /* Use the selection in PB_RT populated by borders to
         * pack the atom position */
        for (int j = 0; j < pb_rt->natoms; j++) {
            int iatom = pb_rt->sel[j];
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

            packbuf_add(pb_r, &r, NULL, NULL);
        }

        if (pb_r->natoms != pb_rt->natoms)
            die("mismatch packed natoms\n");

        dbg("box.%d neigh.%d: packed %d internal ghosts\n",
                box->i, neigh->i, pb_r->natoms);
    }

    box_packbuf_switch(box, PB_R, PB_SEND, PB_PACKING, PB_READY);
    box_packbuf_switch(box, PB_RT, PB_SEND, PB_READING, PB_READY);
}

#pragma oss task label("neigh_border_pack_rt") \
    in(box->r) \
    out({box->pb[PB_RT][PB_SEND][i]->data,    i=0;NNEIGH}) \
    out({box->pb[PB_RT][PB_SEND][i]->natoms,  i=0;NNEIGH}) \
    out({box->pb[PB_RT][PB_SEND][i]->sel,     i=0;NNEIGH})
static void
box_border_pack_rt(Sim *sim, Box *box)
{
    dbg("rank%d.box%d: packing borders with nlocal %d\n",
            sim->rank, box->i, box->nlocal);

    /* Reset all PackBuf from neighbors */
    box_packbuf_clear(box, PB_RT, PB_SEND);
    box_packbuf_switch(box, PB_RT, PB_SEND, PB_READY, PB_PACKING);

    /* Reset ghosts in this box */
    box->nghost = 0;

    int npacked = 0;

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
            PackBuf *pb = box->pb[PB_RT][PB_SEND][neigh->i];

            dbg("rank%d.box%d.atom%d: out of core, sub=%d/%d (%2d %2d %2d), neigh=%d (%2d %2d %2d)\n",
                    sim->rank, box->i, i,
                    j, sub->nneigh,
                    delta[X], delta[Y], delta[Z],
                    neigh->i,
                    neigh->delta[X], neigh->delta[Y], neigh->delta[Z]);

            /* Enforce PBC before packing the atom position */
            if (neigh->wraps) {
                dbg("rank%d.box%d.atom%d: wrapping to neigh %d, position from %e %e %e\n",
                        sim->rank, box->i, i, neigh->i, r[X], r[Y], r[Z]);
    
                for (int d = X; d <= Z; d++)
                    r[d] += neigh->addpbc[d];
    
                dbg("rank%d.box%d.atom%d: wrapped  to neigh %d, position to   %e %e %e\n",
                        sim->rank, box->i, i, neigh->i, r[X], r[Y], r[Z]);
            }

            int type = box->atomtype[i];
            packbuf_add_sel(pb, &r, NULL, &type, i);
        }

        npacked++;
    }

    box_packbuf_switch(box, PB_RT, PB_SEND, PB_PACKING, PB_READY);

    dbg("rank%d.box%d: packed %d atoms\n",
            sim->rank, box->i, npacked);
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
    inout(box->r) \
    out({box->pb[PB_RVT][PB_SEND][i]->data,     i=0;NNEIGH}) \
    out({box->pb[PB_RVT][PB_SEND][i]->natoms,   i=0;NNEIGH})
static void
box_tidy_pack_rvt(Sim *sim, Box *box)
{
    dbg("packing out atoms for box %2d with nlocal %d\n",
            box->i, box->nlocal);

    /* Reset all PackBuf from neighbors */
    box_packbuf_clear(box, PB_RVT, PB_SEND);
    box_packbuf_switch(box, PB_RVT, PB_SEND, PB_READY, PB_PACKING);

    /* Invalidate ghosts, as we are going to modify box->nlocal */
    box->nghost = -666;

    int npacked = 0;

    for (int i = 0; i < box->nlocal; /* nop */) {

        Vec r = { box->r[i][X], box->r[i][Y], box->r[i][Z] };
        int delta[NDIM];

        if (in_domain_delta(r, box->dombox, delta)) {
            i++;
            continue;
        }

        Neigh *neigh = &box->neigh[delta2neigh(delta)];
        PackBuf *pb = &neigh->pb[PB_RVT][PB_SEND];

        /* Enforce PBC before packing the atom position */
        if (neigh->wraps) {
            for (int d = X; d <= Z; d++)
                r[d] += neigh->addpbc[d];
        }

        int type = box->atomtype[i];
        packbuf_add(pb, &r, &box->v[i], &type);
        npacked++;

        /* Fill the hole with one atom from the end */
        int src = box->nlocal - 1, dst = i;
        copy_atom_rvt(box, src, dst);
        box->nlocal--;
    }

    dbg("packed out %d atoms for box %2d with nlocal %d\n",
            npacked, box->i, box->nlocal);

    box_packbuf_switch(box, PB_RVT, PB_SEND, PB_PACKING, PB_READY);
}

void
comm_pack(Sim *sim, enum pb_type type, enum pb_dir dir)
{
    dbg("comm_pack %s\n", PB_TYPENAME(type));

    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];

        switch (type) {
        case PB_R: box_ghost_pack_r(sim, box); break;
        case PB_RT: box_border_pack_rt(sim, box); break;
        case PB_RVT: box_tidy_pack_rvt(sim, box); break;
        default: die("not implemented\n");
        }
    }
}
