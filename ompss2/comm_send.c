#define ENABLE_DEBUG 0
#include "types.h"
#include "log.h"
#include "packbuf.h"

#pragma oss task label("box_send_neigh_task") \
    inout(box->pb[type][neigh->i]->natoms) \
    inout(box->pb[type][neigh->i]->buf)
    //inout(pb->buf) inout(pb->natoms)
static void
box_send_neigh_task(Sim *sim, Box *box, Neigh *neigh, PackBuf *pb,
        enum pb_type type, enum pb_reqtype reqtype)
{
    dbg("%-6s comm_send rank %d, box %d, neigh %d, %s.%s\n",
            "RUN", sim->rank, box->i, neigh->i,
            PB_TYPENAME(type), PB_REQTYPENAME(reqtype));

    packbuf_debug_switch(pb, PB_READY, PB_SENDING);

    packbuf_send(pb, reqtype);

    dbg("%-6s comm_send rank %d, box %d, neigh %d, %s.%s\n",
            "DONE", sim->rank, box->i, neigh->i,
            PB_TYPENAME(type), PB_REQTYPENAME(reqtype));

    packbuf_debug_switch(pb, PB_SENDING, PB_READY);
}

static void
box_send_neigh(Sim *sim, Box *box, Neigh *neigh, enum pb_type type,
        enum pb_reqtype reqtype)
{
    if (neigh->rank == sim->rank) {
        /* No-op: will be copied at recv */
        //dbg("refusing send for rank %d box %d neigh %d\n",
        //        sim->rank, box->i, neigh->i);
        return;
    }

    PackBuf *pb = box->pb[type][neigh->i];

    dbg("%-6s comm_send rank %d, box %d, neigh %d, %s.%s buf=%p buf'=%p\n",
            "CREATE", sim->rank, box->i, neigh->i,
            PB_TYPENAME(type), PB_REQTYPENAME(reqtype),
            &pb->buf, &box->pb[type][neigh->i]->buf);

    box_send_neigh_task(sim, box, neigh, pb, type, reqtype);
}

void
comm_send(Sim *sim, enum pb_type type, enum pb_reqtype reqtype)
{
    //dbg("comm_send %s rank %d\n", PB_TYPENAME(type), sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            box_send_neigh(sim, box, &box->neigh[j], type, reqtype);
        }
    }

    //#pragma oss taskwait
}
