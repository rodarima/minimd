#define ENABLE_DEBUG 0
#include "comm.h"
#include "types.h"
#include "log.h"
#include "packbuf.h"
#include "box.h"
#include "dom.h"
#include "neigh.h"

#include <math.h>

static void
check_max_jump(Sim *sim, Box *box, PackBuf *pb, int i)
{
    if (pb->natoms > pb->nalloc)
        die("%s: natoms=%d larger than nalloc=%d\n",
                pb->name, pb->natoms, pb->nalloc);

    for (int j = 0; j < pb->natoms; j++) {
        int ii = i + j;
        Vec oldr = {
            box->r[ii][X],
            box->r[ii][Y],
            box->r[ii][Z]
        };

        Vec newr = {
            pb->data->buf[j * pb->atomsize + X],
            pb->data->buf[j * pb->atomsize + Y],
            pb->data->buf[j * pb->atomsize + Z]
        };

        /* We should not expect jumps caused by wraps here */

        double dist = sqrt(get_distsq(oldr, newr));
        if (dist > 0.1)
            die("rank %d box %d: jump too large: ii=%d oldr=(%e %e %e) newr=(%e %e %e)\n",
                    sim->rank, box->i, ii,
                    oldr[X], oldr[Y], oldr[Z],
                    newr[X], newr[Y], newr[Z]);

    }
}

static void
neigh_ghost_unpack_r(Sim *sim, Box *box, Neigh *neigh)
{
    PackBuf *pb_r = &neigh->pb[PB_R][PB_RECV];
    PackBuf *pb_rt = &neigh->pb[PB_RT][PB_RECV];

    packbuf_debug_switch(pb_r, PB_READY, PB_RECVING);
    packbuf_debug_switch(pb_rt, PB_READY, PB_READING);

    if (pb_rt->natoms != pb_r->natoms)
        die("ghost_unpack_r: %s mismatch natoms r=%d != rt=%d (neigh natoms=%d)\n",
                pb_r->name, pb_r->natoms, pb_rt->natoms, neigh->recv_natoms_r);

    if (pb_r->natoms != 0) {

        /* Only check the header if we received some atoms */
        if (pb_r->remoterank != sim->rank)
            packbuf_check_header(pb_r);


        int nnew = pb_r->natoms;
        int ntot = box->nlocal + box->nghost + nnew;
        if (ntot > box->nalloc) {
            die("box%d:neigh%d cannot unpack %d atoms, capacity exeeeded\n",
                    box->i, neigh->i, nnew);
        }

        int i = box->nlocal + box->nghost;

        if (ENABLE_MAX_JUMP_CHECK)
            check_max_jump(sim, box, pb_r, i);

        /* The unpack order must be kept the same to match the ghost atom
         * order given by borders */
        packbuf_unpack(pb_r, &box->r[i], NULL, NULL);

        /* We cannot check the domain bounds of the new ghost atom
         * positions, as they are moving around, even exceeding the halo
         * domain */

        box->nghost += pb_r->natoms;
    }

    packbuf_debug_switch(pb_r, PB_RECVING, PB_READY);
    packbuf_debug_switch(pb_rt, PB_READING, PB_READY);
}

/* FIXME: We shouldn't need to use inout */
#pragma oss task label("box_ghost_unpack_r_neigh") \
    inout({box->pb[PB_RT][PB_RECV][i]->natoms,  i=0;NNEIGH}) \
    inout({box->pb[PB_RT][PB_RECV][i]->data,    i=0;NNEIGH}) \
    inout({box->pb[PB_R ][PB_RECV][i]->natoms,  i=0;NNEIGH}) \
    inout({box->pb[PB_R ][PB_RECV][i]->data,    i=0;NNEIGH}) \
    inout(box->r)
static void
box_ghost_unpack_r(Sim *sim, Box *box)
{
    dbg("unpacking r for box %d\n", box->i);

    int old_nghost = box->nghost;
    box->nghost = 0;

    for (int j = 0; j < NNEIGH; j++) {
        /* NOTE: The unpack order must match the neighbor order of the unpack of
         * the rt buffer */
        Neigh *neigh = &box->neigh[j];
        neigh_ghost_unpack_r(sim, box, neigh);
    }

    if (ENABLE_ATOM_COUNT_CHECK) {
        /* Wait until all unpack have finished */
        if (box->nghost != old_nghost) {
            die("nghost atoms don't match %d != %d\n",
                    box->nghost, old_nghost);
        }
    }

    box->fresh_ghost = 0;
}

/* FIXME: We shouldn't need to use inout */
#pragma oss task label("neigh_border_unpack_rt") \
    inout(neigh->pb[PB_RT][PB_RECV].data) \
    inout(neigh->pb[PB_RT][PB_RECV].natoms) \
    inout(neigh->pb[PB_R][PB_RECV].data) \
    inout(neigh->pb[PB_R][PB_RECV].natoms) \
    inout(neigh->recv_natoms_r) \
    out(box->r)
static void
neigh_border_unpack_rt(Sim *sim, Box *box, Neigh *neigh)
{
    PackBuf *pb_rt = &neigh->pb[PB_RT][PB_RECV];
    PackBuf *pb_r = &neigh->pb[PB_R][PB_RECV];

    if (pb_rt->remoterank != sim->rank)
        packbuf_check_header(pb_rt);

    /* Set the natoms to be sent and received in the PB_R buffer */
    neigh->recv_natoms_r = pb_rt->natoms;

    /* Ensure it has room */
    packbuf_grow(pb_r, pb_rt->natoms);
    pb_r->natoms = pb_rt->natoms;

    err("%s setting natoms=%d\n", pb_r->name, pb_r->natoms);

    if (pb_rt->natoms == 0)
        return;

    packbuf_debug_switch(pb_rt, PB_READY, PB_UNPACKING);

    /* Ensure we have room to place the new ghost atoms */
    int nnew = pb_rt->natoms;
    int nend = box->nlocal + box->nghost;
    int ntot = nend + nnew;
    box_realloc(box, ntot);

    /* Unpack the position and type at the end of the local atoms */
    packbuf_unpack(pb_rt, &box->r[nend], NULL, &box->atomtype[nend]);

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

    packbuf_debug_switch(pb_rt, PB_UNPACKING, PB_READY);
}

static void
box_border_unpack_rt(Sim *sim, Box *box)
{
    dbg("unpacking rt for box %d\n", box->i);

    for (int j = 0; j < NNEIGH; j++) {
        Neigh *neigh = &box->neigh[j];
        neigh_border_unpack_rt(sim, box, neigh);
    }

    box->fresh_ghost = 0;
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
    inout(box->r) \
    inout(box->v) \
    inout(box->f) /* May realloc f too */\
    inout(neigh->pb[PB_RVT][PB_RECV].data) \
    inout(neigh->pb[PB_RVT][PB_RECV].natoms)
static void
neigh_tidy_unpack_rvt(Sim *sim, Box *box, Neigh *neigh)
{
    PackBuf *pb = &neigh->pb[PB_RVT][PB_RECV];

    if (pb->remoterank != sim->rank)
        packbuf_check_header(pb);

    packbuf_debug_switch(pb, PB_READY, PB_UNPACKING);

    /* Ensure we have room to place the new local atoms */
    int n = box->nlocal + pb->natoms;
    box_realloc(box, n);

//    /* Here we are receiving new local atoms that have just moved into
//     * our box. Ensure that there was a ghost atom before and that is
//     * not too close to an already existing atom */
//    for (int i = 0; i < pb->natoms; i++) {
//        double *r = &pb->data->buf[i * (NDIM * 2 + 1)];
//        check_atom(sim, box, *(Vec *) r);
//    }

    /* Unpack at the end of the local atoms */
    Vec *r = &box->r[box->nlocal];
    Vec *v = &box->v[box->nlocal];
    int *types = &box->atomtype[box->nlocal];
    packbuf_unpack(pb, r, v, types);

    /* Adjust the number of local atoms in the box */
    box->nlocal = n;
    packbuf_debug_switch(pb, PB_UNPACKING, PB_READY);
}

static void
box_tidy_unpack_rvt(Sim *sim, Box *box)
{
    dbg("unpacking rvt for box %d\n", box->i);

    for (int j = 0; j < NNEIGH; j++) {
        Neigh *neigh = &box->neigh[j];
        neigh_tidy_unpack_rvt(sim, box, neigh);
    }

    box->fresh_ghost = 0;
}

void
comm_unpack(Sim *sim, enum pb_type type, enum pb_dir dir)
{
    dbg("comm_unpack %s\n", PB_TYPENAME(type));

    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];

        switch (type) {
        case PB_R: box_ghost_unpack_r(sim, box); break;
        case PB_RT: box_border_unpack_rt(sim, box); break;
        case PB_RVT: box_tidy_unpack_rvt(sim, box); break;
        default: die("not implemented\n");
        }
    }
}
