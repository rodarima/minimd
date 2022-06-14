#define ENABLE_DEBUG 0
#include "comm.h"
#include "types.h"
#include "log.h"
#include "packbuf.h"
#include "trace.h"

#pragma oss task label("box_send_neigh_task") \
    inout(pb->data) \
    inout(pb->natoms)
static void
box_send_neigh_task(Sim *sim, Box *box, Neigh *neigh, PackBuf *pb,
        enum pb_type type, enum pb_req req)
{
    dbg("%-6s comm_send rank %d, box %d, neigh %d, %s.%s\n",
            "RUN", sim->rank, box->i, neigh->i,
            PB_TYPENAME(type), PB_REQNAME(req));

    char label[1024];
    double t0 = MPI_Wtime();
    sprintf(label, "BEGIN box_send_neigh_task %s transport=%d natoms=%d %s",
        PB_REQNAME(req), pb->transport, pb->natoms, pb->name);
    trace_event(box->i * NNEIGH + neigh->i, label, "#008800");

    if (pb->trace[req].started)
        die("already tracing\n");

    sprintf(pb->trace[req].color, "#ff0000");
    sprintf(pb->trace[req].label, "box_send_neigh_task AFTER WAIT box=%d neigh=%d %s.%s",
            box->i, neigh->i,
            PB_TYPENAME(type), PB_REQNAME(req));

    pb->trace[req].t0 = t0;
    pb->trace[req].h = 0.7;
    pb->trace[req].started = 1;

    packbuf_debug_switch(pb, PB_READY, PB_SENDING);

    /* FIXME: This should be implemented by PackBuf */
    pb->data->header.iter = box->iter;

    if (pb->transport != PB_SHM) {
        packbuf_send(pb, req);
    } else {
        /* No-op: will be copied at recv */
    }

    dbg("%-6s comm_send rank %d, box %d, neigh %d, %s.%s\n",
            "DONE", sim->rank, box->i, neigh->i,
            PB_TYPENAME(type), PB_REQNAME(req));

    packbuf_debug_switch(pb, PB_SENDING, PB_READY);
}

static void
box_send_neigh(Sim *sim, Box *box, Neigh *neigh, enum pb_type type,
        enum pb_dir dir, enum pb_req req)
{
    PackBuf *pb = box->pb[type][dir][neigh->i];

    dbg("%-6s comm_send rank %d, box %d, neigh %d, %s.%s\n",
            "CREATE", sim->rank, box->i, neigh->i,
            PB_TYPENAME(type), PB_REQNAME(req));

    box_send_neigh_task(sim, box, neigh, pb, type, req);
}

void
comm_send(Sim *sim, enum pb_type type, enum pb_dir dir, enum pb_req req)
{
    //dbg("comm_send %s rank %d\n", PB_TYPENAME(type), sim->rank);
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            box_send_neigh(sim, box, &box->neigh[j], type, dir, req);
        }
    }
}
