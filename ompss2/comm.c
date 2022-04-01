#define ENABLE_DEBUG 0
#include "types.h"
#include "log.h"
#include "neigh.h"
#include "dom.h"

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>

#define OFFSETOF(TYPE, ELEMENT) ((size_t)&(((TYPE *)0)->ELEMENT))

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

#pragma oss task label("box_waitmpi_natoms") \
    out({*(char **) &((PackBuf *) (((char *) &box->neigh[i]) + off))->natoms, i=0;NNEIGH}) \
    out({*(char **) &((PackBuf *) (((char *) &box->neigh[i]) + off))->buf,    i=0;NNEIGH})
static void
box_waitmpi_natoms(Sim *sim, Box *box, size_t off, const char *name)
{
    MPI_Request req[NNEIGH];
    int nreq = 0;

    for (int i = 0; i < NNEIGH; i++) {
        Neigh *neigh = &box->neigh[i];
        PackBuf *pb = (PackBuf *) (((char *) neigh) + off);

        packbuf_debug_switch(pb, PB_READY, PB_WAITING);

        if (pb->waitreqn) {
            dbg("rank%d:box%d waiting for natoms in neigh %d\n",
                    sim->rank, box->i, i);
            //MPI_Wait(&pb->reqn, MPI_STATUS_IGNORE);
            memcpy(&req[nreq++], &pb->reqn, sizeof(MPI_Request));
            pb->waitreqn = 0;
        }

    }

    MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE);

    for (int i = 0; i < NNEIGH; i++) {
        Neigh *neigh = &box->neigh[i];
        PackBuf *pb = (PackBuf *) (((char *) neigh) + off);
        packbuf_debug_switch(pb, PB_WAITING, PB_READY);
    }
}

/* FIXME: We should use in() for send buffers */
#pragma oss task label("box_waitmpi_buf") \
    out({*(char **) &((PackBuf *) (((char *) &box->neigh[i]) + off))->natoms, i=0;NNEIGH}) \
    out({*(char **) &((PackBuf *) (((char *) &box->neigh[i]) + off))->buf,    i=0;NNEIGH})
static void
box_waitmpi_buf(Sim *sim, Box *box, size_t off, const char *name)
{
    MPI_Request req[NNEIGH];
    int nreq = 0;

    for (int i = 0; i < NNEIGH; i++) {
        Neigh *neigh = &box->neigh[i];
        PackBuf *pb = (PackBuf *) (((char *) neigh) + off);

        packbuf_debug_switch(pb, PB_READY, PB_WAITING);

        if (pb->waitreq) {
            dbg("rank%d:box%d waiting for buf in neigh %d\n",
                    sim->rank, box->i, i);
            //MPI_Wait(&pb->req, MPI_STATUS_IGNORE);
            memcpy(&req[nreq++], &pb->req, sizeof(MPI_Request));
            pb->waitreq = 0;
        }
    }

    MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE);

    for (int i = 0; i < NNEIGH; i++) {
        Neigh *neigh = &box->neigh[i];
        PackBuf *pb = (PackBuf *) (((char *) neigh) + off);
        packbuf_debug_switch(pb, PB_WAITING, PB_READY);
    }
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
//    dbg("packing out atoms for box %2d with nlocal %d\n",
//            box->i, box->nlocal);

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

//    for (int i = 0; i < NNEIGH; i++) {
//        Neigh *neigh = &box->neigh[i];
//        if (neigh->send_rvt.natoms > 0) {
//            dbg("packed %d atoms in box%d:neigh%02d\n",
//                    neigh->send_rvt.natoms,
//                    box->i, neigh->i);
//        }
//    }
}

static void
box_tidy_send_rvt(Sim *sim, Box *box, Neigh *neigh)
{
    if (neigh->rank != sim->rank) {
        #pragma oss task label("box_tidy_send_rvt:mpisend") \
            in(*(char **)&neigh->send_rvt.buf) \
            in(neigh->send_rvt.natoms)
        {
//            dbg("send_rvt: rank%d:box%d:neigh%d uses mpisend delta=(%d %d %d)\n",
//                    sim->rank, box->i, neigh->i,
//                    neigh->delta[X], neigh->delta[Y], neigh->delta[Z]);

            packbuf_debug_switch(&neigh->send_rvt, PB_READY, PB_SENDING);
            packbuf_mpisend(&neigh->send_rvt);
            packbuf_debug_switch(&neigh->send_rvt, PB_SENDING, PB_READY);
        }
    } else {
//        dbg("send_rvt: rank%d:box%d:neigh%d uses shmcopy delta=(%d %d %d)\n",
//                sim->rank, box->i, neigh->i,
//                neigh->delta[X], neigh->delta[Y], neigh->delta[Z]);

        /* No-op: will be copied in box_tidy_recv_rvt */
    }
}

static void
box_tidy_recv_rvt(Sim *sim, Box *dstbox, Neigh *neigh, int only_natoms)
{
	/* Use MPI for inter process comm */
    if (neigh->rank != sim->rank) {
		if (only_natoms) {
			/* Only natoms */
			#pragma oss task label("box_tidy_recv_rvt:mpirecv:natoms") \
				out(neigh->recv_rvt.recvnatoms) \
				out(*(char **)&neigh->recv_rvt.natoms) \
				out(*(char **)&neigh->recv_rvt.buf) /* Not needed, but
                                                       prevents breaking
                                                       the debug state
                                                     */
			{
				packbuf_debug_switch(&neigh->recv_rvt, PB_READY, PB_RECVING);
				packbuf_mpirecv_natoms(&neigh->recv_rvt);
				packbuf_debug_switch(&neigh->recv_rvt, PB_RECVING, PB_READY);
			}
		} else { /* The buffer */

            /* FIXME: In this task we shouldn't be updating natoms, as
             * it is not yet updated */
			#pragma oss task label("box_tidy_recv_rvt:mpirecv:buf") \
				in(neigh->recv_rvt.recvnatoms) \
				out(*(char **)&neigh->recv_rvt.buf) \
				out(*(char **)&neigh->recv_rvt.natoms)
			{
				packbuf_debug_switch(&neigh->recv_rvt, PB_READY, PB_RECVING);
				packbuf_mpirecv_buf(&neigh->recv_rvt, neigh->recv_rvt.recvnatoms);
				packbuf_debug_switch(&neigh->recv_rvt, PB_RECVING, PB_READY);
			}
		}
    } else if (!only_natoms) { /* Only shmcopy once, in the last call */

        /* Shared memory for intra-process. This can be avoided if
         * we pack directly into the receiving buffer. */
        Neigh *dstneigh = neigh;
        Box *srcbox = dstneigh->box;
        Neigh *srcneigh = &srcbox->neigh[opposite_neigh(dstneigh->i)];

        #pragma oss task label("box_tidy_recv_rvt:shmcopy") \
            in(*(char **)&srcneigh->send_rvt.buf) \
            in(srcneigh->send_rvt.natoms) \
            out(*(char **)&dstneigh->recv_rvt.buf) \
            out(dstneigh->recv_rvt.natoms)
        {
            packbuf_debug_switch(&srcneigh->send_rvt, PB_READY, PB_COPYING);
            packbuf_debug_switch(&dstneigh->recv_rvt, PB_READY, PB_COPYING);

//            dbg("rank%d:box%d:neigh%02d shmcopy\n",
//                    sim->rank, dstbox->i, dstneigh->i);
            packbuf_shmcopy(&srcneigh->send_rvt, &dstneigh->recv_rvt);

            packbuf_debug_switch(&srcneigh->send_rvt, PB_COPYING, PB_READY);
            packbuf_debug_switch(&dstneigh->recv_rvt, PB_COPYING, PB_READY);
        }

//        if (srcneigh->send_rvt.natoms > 0) {
//            dbg("shmcopy %d atoms from box%d:neigh%02d -> box%d:neigh%02d\n",
//                    srcneigh->send_rvt.natoms,
//                    srcbox->i, srcneigh->i,
//                    dstbox->i, dstneigh->i);
//        }
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

#pragma oss task label("box_tidy_unpack_rvt") \
    inout(*(char **)&box->r) \
    inout(*(char **)&box->v) \
    inout(*(char **)&box->f) /* May realloc f too */\
    inout(*(char **)&neigh->recv_rvt.buf) \
    inout(*(char **)&neigh->recv_rvt.natoms)
static void
box_tidy_unpack_rvt(Sim *sim, Box *box, Neigh *neigh)
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

void
comm_tidy(Sim *sim)
{
    dbg("rank%d -- comm_tidy -- begins\n", sim->rank);

    if (ENABLE_NONBLOCKING_MPI) {
        dbg("rank%d -- comm_tidy -- waitmpi send_rvt buffer\n", sim->rank);
        for (int i = 0; i < sim->nboxes; i++) {
            box_waitmpi_natoms(sim, &sim->box[i], OFFSETOF(Neigh, send_rvt),
				"comm_tidy -- waitmpi send_rvt natoms");
            box_waitmpi_buf(sim, &sim->box[i], OFFSETOF(Neigh, send_rvt),
				"comm_tidy -- waitmpi send_rvt buf");
        }
    }

    dbg("rank%d -- comm_tidy -- pack\n", sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        box_tidy_pack_rvt(sim, box);
    }

    dbg("rank%d -- comm_tidy -- send send_rvt natoms+buf\n", sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int i = 0; i < NNEIGH; i++) {
            Neigh *neigh = &box->neigh[i];
            box_tidy_send_rvt(sim, box, neigh);
        }
    }

    dbg("rank%d -- comm_tidy -- recv recv_rvt natoms\n", sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int i = 0; i < NNEIGH; i++) {
            Neigh *neigh = &box->neigh[i];
            box_tidy_recv_rvt(sim, box, neigh->opposite, 1);
        }
    }

    /* Wait for all natom messages to get buffer sizes */
    if (ENABLE_NONBLOCKING_MPI) {
        dbg("rank%d -- comm_tidy -- waitmpi recv_rvt natoms\n", sim->rank);
        for (int i = 0; i < sim->nboxes; i++) {
            box_waitmpi_natoms(sim, &sim->box[i], OFFSETOF(Neigh, recv_rvt),
				"comm_tidy -- waitmpi recv_rvt natoms");
		}
    }

    dbg("rank%d -- comm_tidy -- recv recv_rvt buf\n", sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int i = 0; i < NNEIGH; i++) {
            Neigh *neigh = &box->neigh[i];
            box_tidy_recv_rvt(sim, box, neigh->opposite, 0);
        }
    }

    if (ENABLE_NONBLOCKING_MPI) {
        dbg("rank%d -- comm_tidy -- waitmpi recv_rvt buf\n", sim->rank);
        for (int i = 0; i < sim->nboxes; i++) {
            box_waitmpi_buf(sim, &sim->box[i], OFFSETOF(Neigh, recv_rvt),
				"comm_tidy -- waitmpi recv_rvt buf");
		}
    }

    dbg("rank%d -- comm_tidy -- unpack recv_rvt\n", sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int i = 0; i < NNEIGH; i++) {
            Neigh *neigh = &box->neigh[i];
            box_tidy_unpack_rvt(sim, box, neigh->opposite);
        }
    }

    dbg("rank%d -- comm_tidy -- ends\n", sim->rank);
}

#pragma oss task label("box_border_pack_rt") \
    in(*(char **)&box->r) \
	out({box->neigh[i].send_rt.buf, i=0;NNEIGH}) \
	out({box->neigh[i].send_rt.natoms, i=0;NNEIGH})
static void
box_border_pack_rt(Sim *sim, Box *box)
{
//    dbg("packing borders for box %2d with nlocal %d\n",
//            box->i, box->nlocal);

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

//            dbg("atom %d out of core, delta sub (%2d %2d %2d), neigh %d/%d\n",
//                    i, delta[X], delta[Y], delta[Z], neigh->i,
//                    sub->nneigh);

            /* Enforce PBC before packing the atom position */
            if (neigh->wraps) {
//                dbg("wrapping neigh %d atom %d position from %e %e %e\n",
//                        neigh->i, i, r[X], r[Y], r[Z]);
    
                for (int d = X; d <= Z; d++)
                    r[d] += neigh->addpbc[d];
    
//                dbg("wrapped  neigh %d atom %d position to   %e %e %e\n",
//                        neigh->i, i, r[X], r[Y], r[Z]);
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
}

static void
box_border_send_rt(Sim *sim, Box *box, Neigh *neigh)
{
    if (neigh->rank != sim->rank) {
        #pragma oss task label("box_border_send_rt:mpisend") \
            in(*(char **)&neigh->send_rt.buf) \
            in(neigh->send_rt.natoms)
        {
            packbuf_debug_switch(&neigh->send_rt, PB_READY, PB_SENDING);
            packbuf_mpisend(&neigh->send_rt);
            packbuf_debug_switch(&neigh->send_rt, PB_SENDING, PB_READY);
        }
    } else {
        /* No-op: will be copied via shared memory at recv */
    }
}

static void
box_border_recv_rt(Sim *sim, Box *box, Neigh *dstneigh, int only_natoms)
{
    if (dstneigh->rank != sim->rank) {
		if (only_natoms) {
			#pragma oss task label("box_border_recv_rt:mpirecv:natoms") \
				out(dstneigh->recv_rt.recvnatoms) \
				out(*(char **)&dstneigh->recv_rt.buf) \
				out(*(char **)&dstneigh->recv_rt.natoms)
			packbuf_mpirecv_natoms(&dstneigh->recv_rt);
		} else {
			#pragma oss task label("box_border_recv_rt:mpirecv:buf") \
				in(dstneigh->recv_rt.recvnatoms) \
				out(*(char **)&dstneigh->recv_rt.buf) \
				out(*(char **)&dstneigh->recv_rt.natoms)
			packbuf_mpirecv_buf(&dstneigh->recv_rt, dstneigh->recv_rt.recvnatoms);
		}
    } else if (!only_natoms) {
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
        #pragma oss task label("box_border_recv_rt:shmcopy") \
            in(*(char **)&srcneigh->send_rt.buf) \
            in(srcneigh->send_rt.natoms) \
			out(*(char **)&dstneigh->recv_rt.buf) \
			out(dstneigh->recv_rt.natoms)
        {
            packbuf_debug_switch(&srcneigh->send_rt, PB_READY, PB_COPYING);
            packbuf_debug_switch(&dstneigh->recv_rt, PB_READY, PB_COPYING);
            /* Clear receive buffer */
            //dbg("clearing box%d:neigh%d recv_rt\n", box->i, dstneigh->i);
            packbuf_clear(&dstneigh->recv_rt);
            packbuf_shmcopy(&srcneigh->send_rt, &dstneigh->recv_rt);
            //dbg("set box%d:neigh%d recv_rt %d\n",
            //        box->i, dstneigh->i, dstneigh->recv_rt.natoms);
            //dbg("box %d neigh %d has in recv_rt %d atoms\n",
            //        box->i, dstneigh->i, dstneigh->recv_rt.natoms);
            packbuf_debug_switch(&dstneigh->recv_rt, PB_COPYING, PB_READY);
            packbuf_debug_switch(&srcneigh->send_rt, PB_COPYING, PB_READY);
        }
    }
}

/* FIXME: We shouldn't need to use inout */
#pragma oss task label("box_border_unpack_rt") \
    inout(*(char **)&neigh->recv_rt.buf) \
    inout(*(char **)&neigh->recv_rt.natoms) \
	out(*(char **)&box->r)
static void
box_border_unpack_rt(Sim *sim, Box *box, Neigh *neigh)
{
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

void
comm_borders(Sim *sim)
{
    if (ENABLE_NONBLOCKING_MPI) {
        dbg("rank%d -- comm_borders -- waitmpi send_rt buf+natoms\n", sim->rank);
        for (int i = 0; i < sim->nboxes; i++) {
            box_waitmpi_natoms(sim, &sim->box[i], OFFSETOF(Neigh, send_rt),
					"comm_borders -- waitmpi send_rt natoms");
			box_waitmpi_buf(sim, &sim->box[i], OFFSETOF(Neigh, send_rt),
					"comm_borders -- waitmpi send_rt buf");
        }
    }

    dbg("rank%d -- comm_borders -- pack send_rt buf+natoms\n", sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        box_border_pack_rt(sim, box);
    }

    dbg("rank%d -- comm_borders -- send send_rt buf+natoms\n", sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            Neigh *neigh = &box->neigh[j];
            box_border_send_rt(sim, box, neigh);
        }
    }

    dbg("rank%d -- comm_borders -- recv recv_rt natoms\n", sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            /* Use the inverse order for reception */
            Neigh *dstneigh = box->neigh[j].opposite;
            box_border_recv_rt(sim, box, dstneigh, 1);
        }
    }

    if (ENABLE_NONBLOCKING_MPI) {
        dbg("rank%d -- comm_borders -- waitmpi recv_rt natoms\n", sim->rank);
        for (int i = 0; i < sim->nboxes; i++) {
            box_waitmpi_natoms(sim, &sim->box[i], OFFSETOF(Neigh, recv_rt),
					"comm_borders -- waitmpi recv_rt natoms");
		}
    }

    dbg("rank%d -- comm_borders -- recv recv_rt buf\n", sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            /* Use the inverse order for reception */
            Neigh *dstneigh = box->neigh[j].opposite;
            box_border_recv_rt(sim, box, dstneigh, 0);
        }
    }

    if (ENABLE_NONBLOCKING_MPI) {
        dbg("rank%d -- comm_borders -- waitmpi recv_rt buf\n", sim->rank);
        for (int i = 0; i < sim->nboxes; i++) {
            box_waitmpi_buf(sim, &sim->box[i], OFFSETOF(Neigh, recv_rt),
					"comm_borders -- waitmpi recv_rt buf");
		}
    }

    dbg("rank%d -- comm_borders -- unpack\n", sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            /* Use the inverse order for unpack */
            Neigh *opp = box->neigh[j].opposite;
            box_border_unpack_rt(sim, box, opp);
        }
        //if (box->nghost == 0)
        //    abort();
    }

    dbg("rank%d -- comm_borders -- ends\n", sim->rank);
}

#pragma oss task label("box_ghost_pack_r") \
    in(*(char **)&box->r) \
    in({*(char **)&box->neigh[i].send_rt.buf, i=0;NNEIGH}) \
    inout({*(char **)&box->neigh[i].send_r.buf, i=0;NNEIGH})
static void
box_ghost_pack_r(Sim *sim, Box *box)
{
//    dbg("packing internal ghosts from box %d\n", box->i);

    /* Reset all PackBuf from neighbors */
    for (int i = 0; i < NNEIGH; i++) {
        packbuf_clear(&box->neigh[i].send_r);
        packbuf_debug_switch(&box->neigh[i].send_r, PB_READY, PB_PACKING);
    }

    /* Use the selection in send_rt populated by borders to pack the
     * atom position */
    for (int i = 0; i < NNEIGH; i++) {
        Neigh *neigh = &box->neigh[i];

        //packbuf_debug_switch(&neigh->send_rt, PB_READY, PB_READING);
        PackBuf *pb = &neigh->send_rt;

        for (int j = 0; j < pb->natoms; j++) {
            int iatom = pb->sel[j];
            if (iatom < 0 || iatom >= box->nlocal) {
                die("atom %d outside local range\n", iatom);
            }

            Vec r = { box->r[iatom][X], box->r[iatom][Y], box->r[iatom][Z] };

            /* Enforce PBC before packing the atom position */
            if (neigh->wraps) {
                for (int d = X; d <= Z; d++)
                    r[d] += neigh->addpbc[d];
            }

            packbuf_add(&neigh->send_r, &r, NULL, NULL);
        }

        //dbg("box %d neigh %d: packed %d internal ghosts\n",
        //        box->i, neigh->i, neigh->send_r.natoms);

        //packbuf_debug_switch(&neigh->send_rt, PB_READING, PB_READY);
        packbuf_debug_switch(&neigh->send_r, PB_PACKING, PB_READY);
    }
}

static void
box_ghost_send_r(Sim *sim, Box *box, Neigh *srcneigh)
{
    if (srcneigh->rank != sim->rank) {
        /* Only send the atom buffer, the receive end already
         * knows the size */
        #pragma oss task label("box_ghost_send_r") \
            in(*(char **)&srcneigh->send_r.buf) \
            in(srcneigh->send_r.natoms)
        {
            packbuf_debug_switch(&srcneigh->send_r, PB_READY, PB_PACKING);
            packbuf_mpisend_buf(&srcneigh->send_r);
            packbuf_debug_switch(&srcneigh->send_r, PB_PACKING, PB_READY);
        }
    } else {
        /* No-op: will be copied via shared memory at recv */
    }
}

static void
box_ghost_recv_r(Sim *sim, Box *box, Neigh *dstneigh)
{
    if (dstneigh->rank != sim->rank) {
        #pragma oss task label("box_ghost_recv_r:mpirecv") \
            firstprivate(dstneigh) \
            in(*(char **)&dstneigh->recv_rt.buf) /* For natoms only */ \
            in(*(char **)&dstneigh->recv_rt.natoms) \
            out(*(char **)&dstneigh->recv_r.buf) \
            out(*(char **)&dstneigh->recv_r.natoms)
        {
            //packbuf_debug_switch(&dstneigh->recv_rt, PB_READY, PB_READING);
            packbuf_debug_switch(&dstneigh->recv_r, PB_READY, PB_RECVING);
            /* Get the number of atoms to be received from the pack
             * buffer used in the borders */
            dbg("reading box%d:neigh%d recv_rt\n", box->i, dstneigh->i);
            int natoms = dstneigh->recv_rt.natoms;
            packbuf_clear(&dstneigh->recv_r);
            packbuf_mpirecv_buf(&dstneigh->recv_r, natoms);

            //packbuf_debug_switch(&dstneigh->recv_rt, PB_READING, PB_READY);
            packbuf_debug_switch(&dstneigh->recv_r, PB_RECVING, PB_READY);
        }
    } else {
        /* Get the source box from the neighbor and find the
         * neighbor which contains the send buffer */
        Box *srcbox = dstneigh->box;
        Neigh *srcneigh = &srcbox->neigh[opposite_neigh(dstneigh->i)];

        #pragma oss task label("box_ghost_recv_r:shmcopy") \
            firstprivate(srcneigh, dstneigh) \
            in(*(char **)&dstneigh->recv_rt.buf) /* For natoms only */ \
            in(*(char **)&srcneigh->send_r.buf) \
            out(*(char **)&dstneigh->recv_r.buf)
        {
            //packbuf_debug_switch(&dstneigh->recv_rt, PB_READY, PB_READING);
            packbuf_debug_switch(&dstneigh->recv_r,  PB_READY, PB_COPYING);
            packbuf_debug_switch(&srcneigh->send_r,  PB_READY, PB_COPYING);

            //dbg("reading box%d:neigh%d recv_rt\n", box->i, dstneigh->i);
            if (srcneigh->send_r.natoms != dstneigh->recv_rt.natoms) {
                dbg("box %d srcneigh %d dstneigh %d: natoms don't match\n"
                            "  srcneigh->send_r.natoms = %d != dstneigh->recv_rt.natoms = %d\n",
                        box->i, srcneigh->i, dstneigh->i,
                        srcneigh->send_r.natoms,
                        dstneigh->recv_rt.natoms);
                sleep(1);
                abort();
            }
            packbuf_clear(&dstneigh->recv_r);
            packbuf_shmcopy(&srcneigh->send_r, &dstneigh->recv_r);
            //dbg("setting box%d:neigh%d recv_r natoms=%d\n",
            //        box->i, dstneigh->i, dstneigh->recv_r.natoms);

            //packbuf_debug_switch(&dstneigh->recv_rt, PB_READING, PB_READY);
            packbuf_debug_switch(&dstneigh->recv_r,  PB_COPYING, PB_READY);
            packbuf_debug_switch(&srcneigh->send_r,  PB_COPYING, PB_READY);
        }
    }

//  if (natoms != 0) {
//      dbg("box %d: unpacked %d ghost atoms from neigh %d\n",
//              box->i, natoms, j);
//  }
}

static void
box_ghost_unpack_r_neigh(Sim *sim, Box *box, Neigh *neigh)
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
    in({*(char **)&box->neigh[i].recv_rt.natoms,    i=0;NNEIGH}) \
    in({*(char **)&box->neigh[i].recv_rt.buf,       i=0;NNEIGH}) \
    in({*(char **)&box->neigh[i].recv_r.natoms,     i=0;NNEIGH}) \
    in({*(char **)&box->neigh[i].recv_r.buf,        i=0;NNEIGH}) \
    inout(*(char **)&box->r)
static void
box_ghost_unpack_r(Sim *sim, Box *box)
{
    int old_nghost = box->nghost;
    box->nghost = 0;

    for (int j = 0; j < NNEIGH; j++) {
        /* Use the inverse order for unpack */
        Neigh *opp = box->neigh[j].opposite;
        box_ghost_unpack_r_neigh(sim, box, opp);
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

/* Send/recv ghost positions (communicate) */
void
comm_ghost_position(Sim *sim)
{
    /* No need to wait for natoms as its not sent */
    if (ENABLE_NONBLOCKING_MPI) {
        dbg("rank%d -- comm_ghost -- waitmpi send_r buf\n", sim->rank);
        for (int i = 0; i < sim->nboxes; i++) {
            box_waitmpi_buf(sim, &sim->box[i], OFFSETOF(Neigh, send_r),
					"comm_ghost -- waitmpi send_r buf");
		}
    }

    dbg("rank%d -- comm_ghost -- pack_r buf\n", sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        box_ghost_pack_r(sim, box);
    }

    dbg("rank%d -- comm_ghost -- send_r buf\n", sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            Neigh *srcneigh = &box->neigh[j];
            box_ghost_send_r(sim, box, srcneigh);
        }
    }

    /* No need to receive the number of atoms, as it is known */
    dbg("rank%d -- comm_ghost -- recv_r buf\n", sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            /* Use the inverse order for reception */
            Neigh *dstneigh = box->neigh[j].opposite;
            box_ghost_recv_r(sim, box, dstneigh);
        }
    }

    /* Wait for all buffer messages to complete before unpack */
    if (ENABLE_NONBLOCKING_MPI) {
        dbg("rank%d -- comm_ghost -- waitmpi recv_r buf\n", sim->rank);
        for (int i = 0; i < sim->nboxes; i++) {
            box_waitmpi_buf(sim, &sim->box[i], OFFSETOF(Neigh, recv_r),
					"comm_ghost -- waitmpi recv_r buf");
		}
    }


    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        box_ghost_unpack_r(sim, box);
    }
}
