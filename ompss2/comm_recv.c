#define ENABLE_DEBUG 0
#include "types.h"
#include "log.h"
#include "packbuf.h"
#include "neigh.h"

/** Exchanges the PackBuf of the given type and req by using a internode
 * communication task. */
static void
neigh_recv_internode(Sim *sim, Box *box, Neigh *neigh,
        enum pb_type type, enum pb_req reqtype)
{
    PackBuf *pb = box->pb[neigh->i][type][PB_RECV];

    dbg("%-6s comm_recv rank %d, box %d, neigh %d, %s.%s\n",
            "CREATE", sim->rank, box->i, neigh->i,
            PB_TYPENAME(type), PB_REQNAME(reqtype));

    /* TODO: We should move the task dependencies to the functions that modify
     * the actual data, rather than here */

    #pragma oss task label("neigh_recv_internode") \
        inout(pb->data) inout(pb->natoms)
    {
        dbg("%-6s comm_recv rank %d, box %d, neigh %d, %s.%s\n",
                "RUN", sim->rank, box->i, neigh->i,
                PB_TYPENAME(type), PB_REQNAME(reqtype));

        packbuf_debug_switch(pb, PB_READY, PB_RECVING);
        packbuf_recv(pb, reqtype);

        dbg("%-6s comm_recv rank %d, box %d, neigh %d, %s.%s\n",
                "DONE", sim->rank, box->i, neigh->i,
                PB_TYPENAME(type), PB_REQNAME(reqtype));

        packbuf_debug_switch(pb, PB_RECVING, PB_READY);
    }
}

static void
neigh_recv_intranode(Sim *sim, Box *box, Neigh *neigh,
        enum pb_type type, enum pb_req reqtype)
{
    if (reqtype == PB_NATOMS) {
        /* No-op as we already know the size */
        //dbg("noop for box %d neigh %d with reqtype=%s\n",
        //        box->i, neigh->i, PB_REQNAME(reqtype));
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

    PackBuf *send_pb = send_box->pb[send_idir][type][PB_SEND];
    PackBuf *recv_pb = recv_box->pb[recv_idir][type][PB_RECV];

    #pragma oss task label("neigh_recv_intranode:shmcopy") \
        inout(send_pb->data) \
        inout(send_pb->natoms) \
        inout(recv_pb->data) \
        inout(recv_pb->natoms)
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
        enum pb_type type, enum pb_req reqtype)
{
    /* Use MPI for inter process comm */
    if (neigh->rank != sim->rank) {
        neigh_recv_internode(sim, box, neigh, type, reqtype);
    } else {
        neigh_recv_intranode(sim, box, neigh, type, reqtype);
    }
}

void
comm_recv(Sim *sim, enum pb_type type, enum pb_req reqtype)
{
    dbg("comm_recv type=%s reqtype=%s\n",
            PB_TYPENAME(type), PB_REQNAME(reqtype));

    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++)
            neigh_recv(sim, box, &box->neigh[j], type, reqtype);
    }
}
