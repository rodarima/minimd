#define ENABLE_DEBUG 0
#include "types.h"
#include "log.h"
#include "packbuf.h"
#include "neigh.h"

static void
neigh_recv_internode(Sim *sim, Box *box, Neigh *neigh,
        enum pb_type type, enum pb_reqtype reqtype)
{
    PackBuf *pb = box->pb[type][neigh->i];

    /* TODO: We should move the task dependencies to the functions that modify
     * the actual data, rather than here */

    if (reqtype == PB_NATOMS) {
        #pragma oss task label("neigh_recv_internode:natoms") \
            out(pb->recvnatoms)
        {
            packbuf_debug_switch(pb, PB_READY, PB_RECVING);
            packbuf_recv(pb, PB_NATOMS);
            packbuf_debug_switch(pb, PB_RECVING, PB_READY);
        }
    } else {
        if (ENABLE_NONBLOCKING_MPI) {
            #pragma oss task label("neigh_recv_internode:nb:buf") \
                in(pb->recvnatoms) \
                out(pb->buf) /* pb->natoms will be updated in the MPI_Wait */
            {
                packbuf_debug_switch(pb, PB_READY, PB_RECVING);
                packbuf_recv(pb, PB_BUF);
                packbuf_debug_switch(pb, PB_RECVING, PB_READY);
            }
        } else {
            #pragma oss task label("neigh_recv_internode:b:buf") \
                in(pb->recvnatoms) \
                out(pb->natoms) \
                out(pb->buf)
            {
                packbuf_debug_switch(pb, PB_READY, PB_RECVING);
                packbuf_recv(pb, PB_BUF);
                packbuf_debug_switch(pb, PB_RECVING, PB_READY);
            }
        }
    }

}

static void
neigh_recv_intranode(Sim *sim, Box *box, Neigh *neigh,
        enum pb_type type, enum pb_reqtype reqtype)
{
    if (reqtype == PB_NATOMS) {
        /* No-op as we already know the size */
        //dbg("noop for box %d neigh %d with reqtype=%s\n",
        //        box->i, neigh->i, PB_REQTYPENAME(reqtype));
        return;
    }

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

    int recv_idir = neigh->i;
    int send_idir = opposite_neigh(recv_idir);

    Box *recv_box = box;
    Box *send_box = neigh->box;

    enum pb_type sendtype = packbuf_opposite_dir(type);

    PackBuf *send_pb = send_box->pb[sendtype][send_idir];
    PackBuf *recv_pb = recv_box->pb[type][recv_idir];

    #pragma oss task label("neigh_recv_intranode:shmcopy") \
        in(send_pb->buf) in(send_pb->natoms) \
        out(recv_pb->buf) out(recv_pb->natoms)
    {
        //dbg("shmcopy box%d:neigh%d:pb.%p --(%d)--> box%d:neigh%d:pb.%p\n",
        //        send_box->i, send_idir,
        //        send_pb,
        //        send_pb->natoms,
        //        recv_box->i, recv_idir, recv_pb);

        packbuf_debug_switch(send_pb, PB_READY, PB_COPYING);
        packbuf_debug_switch(recv_pb, PB_READY, PB_COPYING);

        /* Clear receive buffer */
        packbuf_clear(recv_pb);
        packbuf_shmcopy(send_pb, recv_pb, reqtype);

        packbuf_debug_switch(recv_pb, PB_COPYING, PB_READY);
        packbuf_debug_switch(send_pb, PB_COPYING, PB_READY);
    }
}

static void
neigh_recv(Sim *sim, Box *box, Neigh *neigh,
        enum pb_type type, enum pb_reqtype reqtype)
{
    /* Use MPI for inter process comm */
    if (neigh->rank != sim->rank) {
        neigh_recv_internode(sim, box, neigh, type, reqtype);
    } else {
        neigh_recv_intranode(sim, box, neigh, type, reqtype);
    }
}

void
comm_recv(Sim *sim, enum pb_type type, enum pb_reqtype reqtype)
{
    dbg("comm_recv type=%s reqtype=%s\n",
            PB_TYPENAME(type), PB_REQTYPENAME(reqtype));

    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++)
            neigh_recv(sim, box, &box->neigh[j], type, reqtype);
    }
}
