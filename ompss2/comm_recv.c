#define ENABLE_DEBUG 0
#include "comm.h"
#include "types.h"
#include "log.h"
#include "packbuf.h"
#include "neigh.h"

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
//    if (req == PB_NATOMS) {
//        /* No-op as we already know the size */
//        return;
//    }

    /* FIXME: This is a violation of concerns. The PackBuf module should be the
     * only one which knows which other PB are being used and place the
     * dependencies accordingly. */

    PackBuf *recv_pb = box->pb[type][PB_RECV][neigh->i];
    PackBuf *send_pb = recv_pb->shm.remote;

    #pragma oss task label("neigh_recv_intranode:shmcopy") \
        inout(send_pb->data) \
        inout(send_pb->natoms) \
        inout(recv_pb->data) \
        inout(recv_pb->natoms)
    {
        packbuf_debug_switch(send_pb, PB_READY, PB_COPYING);
        packbuf_debug_switch(recv_pb, PB_READY, PB_COPYING);

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
