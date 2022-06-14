#define ENABLE_DEBUG 0
#include "comm.h"
#include "types.h"
#include "log.h"
#include "packbuf.h"
#include "neigh.h"
#include "trace.h"

/** Exchanges the PackBuf of the given type and req by using a internode
 * communication task. */
static void
neigh_recv_internode(Sim *sim, Box *box, Neigh *neigh,
        enum pb_type type, enum pb_req req)
{
    PackBuf *pb = box->pb[type][PB_RECV][neigh->i];

    dbg("%-6s comm_recv rank %d, box %d, neigh %d, %s.%s\n",
            "CREATE", sim->rank, box->i, neigh->i,
            PB_TYPENAME(type), PB_REQNAME(req));

    /* TODO: We should move the task dependencies to the functions that modify
     * the actual data, rather than here */

    #pragma oss task label("neigh_recv_internode") \
        inout(pb->data) inout(pb->natoms)
    {
        dbg("%-6s comm_recv rank %d, box %d, neigh %d, %s.%s\n",
                "RUN", sim->rank, box->i, neigh->i,
                PB_TYPENAME(type), PB_REQNAME(req));

        char label[1024];
        double t0 = MPI_Wtime();
        sprintf(label, "BEGINS neigh_recv_internode %s natoms=%d %s",
            PB_REQNAME(req), pb->natoms, pb->name);
        trace_event(box->i * NNEIGH + neigh->i, label, "#008800");

        if (pb->trace[req].started)
            die("already tracing\n");

        sprintf(pb->trace[req].color, "#0000ff");
        sprintf(pb->trace[req].label, "ENDS comm_recv internode box=%d neigh=%d %s.%s",
                box->i, neigh->i,
                PB_TYPENAME(type), PB_REQNAME(req));

        pb->trace[req].t0 = MPI_Wtime();
        pb->trace[req].h = 0.5;
        pb->trace[req].started = 1;

        packbuf_debug_switch(pb, PB_READY, PB_RECVING);
        packbuf_recv(pb, req);

        dbg("%-6s comm_recv rank %d, box %d, neigh %d, %s.%s\n",
                "DONE", sim->rank, box->i, neigh->i,
                PB_TYPENAME(type), PB_REQNAME(req));

        packbuf_debug_switch(pb, PB_RECVING, PB_READY);
    }
}

static void
neigh_recv_intranode(Sim *sim, Box *box, Neigh *neigh,
        enum pb_type type, enum pb_req req)
{
    if (req == PB_NATOMS) {
        /* No-op as we already know the size */
        //dbg("noop for box %d neigh %d with req=%s\n",
        //        box->i, neigh->i, PB_REQNAME(req));
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

    /* FIXME: This is a violation of concerns. The PackBuf module should be the
     * only one which knows which other PB are being used and place the
     * dependencies accordingly. */

    PackBuf *send_pb = send_box->pb[type][PB_SEND][send_idir];
    PackBuf *recv_pb = recv_box->pb[type][PB_RECV][recv_idir];

    if (recv_pb->shm.remote == NULL)
        die("%s: remote pointer is NULL\n",
                recv_pb->name);

    if (send_pb->shm.remote != NULL)
        die("%s: send remote pointer is not NULL\n",
                send_pb->name);

    if (recv_pb->shm.remote != send_pb)
        die("remote pointer mismatch: recv remote %s != send %s\n",
                recv_pb->shm.remote->name,
                send_pb->name);

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

        PackBuf *pb = recv_pb;

        char label[1024];
        double t0 = MPI_Wtime();
        sprintf(label, "BEGINS neigh_recv_intranode %s %s", PB_REQNAME(req), pb->name);
        trace_event(box->i * NNEIGH + neigh->i, label, "#0088ff");

        packbuf_debug_switch(send_pb, PB_READY, PB_COPYING);
        packbuf_debug_switch(recv_pb, PB_READY, PB_COPYING);

        if (pb->trace[req].started)
            die("already tracing\n");

        sprintf(pb->trace[req].color, "#00ffff");
        sprintf(pb->trace[req].label, "ENDS comm_recv intranode box=%d neigh=%d %s.%s",
                box->i, neigh->i,
                PB_TYPENAME(type), PB_REQNAME(req));

        pb->trace[req].t0 = t0;
        pb->trace[req].h = 0.3;
        pb->trace[req].started = 1;

        /* Clear receive buffer */
        packbuf_clear(recv_pb);
        packbuf_recv(recv_pb, req);

        packbuf_debug_switch(recv_pb, PB_COPYING, PB_READY);
        packbuf_debug_switch(send_pb, PB_COPYING, PB_READY);
    }
}

static void
neigh_recv(Sim *sim, Box *box, Neigh *neigh,
        enum pb_type type, enum pb_dir dir, enum pb_req req)
{

    /* Use MPI for inter process comm */
    if (neigh->rank != sim->rank) {
        neigh_recv_internode(sim, box, neigh, type, req);
    } else {
        neigh_recv_intranode(sim, box, neigh, type, req);
    }
}

void
comm_recv(Sim *sim,
        enum pb_type type, enum pb_dir dir, enum pb_req req)
{
    dbg("comm_recv type=%s req=%s\n",
            PB_TYPENAME(type), PB_REQNAME(req));

    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++)
            neigh_recv(sim, box, &box->neigh[j], type, dir, req);
    }
}
