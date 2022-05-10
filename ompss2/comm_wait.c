#define ENABLE_DEBUG 0
#include "comm.h"
#include "types.h"
#include "log.h"
#include "packbuf.h"

static void
box_waitmpi_task(Sim *sim, Box *box,
        enum pb_type type, enum pb_dir dir, enum pb_req req)
{
    dbg("%-6s comm_wait rank %d, box %d, %s.%s\n",
            "RUN", sim->rank, box->i,
            PB_TYPENAME(type), PB_REQNAME(req));

    for (int i = 0; i < NNEIGH; i++) {
        PackBuf *pb = box->pb[type][dir][i];
        packbuf_debug_switch(pb, PB_READY, PB_WAITING);
    }

    packbuf_mpi_waitn(box->pb[type][dir], NNEIGH, req);

    dbg("%-6s comm_wait rank %d, box %d, %s.%s\n",
            "DONE", sim->rank, box->i,
            PB_TYPENAME(type), PB_REQNAME(req));

    for (int i = 0; i < NNEIGH; i++) {
        PackBuf *pb = box->pb[type][dir][i];
        packbuf_debug_switch(pb, PB_WAITING, PB_READY);
    }
}

static void
box_waitmpi(Sim *sim, Box *box,
        enum pb_type type, enum pb_dir dir, enum pb_req req)
{

    #pragma oss task label("box_waitmpi") \
        inout({box->pb[type][dir][i]->natoms, i=0;NNEIGH}) \
        inout({box->pb[type][dir][i]->data,   i=0;NNEIGH})
    {
        box_waitmpi_task(sim, box, type, dir, req);
    }
}

void
comm_wait(Sim *sim,
        enum pb_type type, enum pb_dir dir, enum pb_req req)
{
    if (!NEED_EXPLICIT_WAIT)
        return;

    //dbg("rank%d -- comm_tidy -- recv recv_rvt buf\n", sim->rank);

    for (int i = 0; i < sim->nboxes; i++) {
        dbg("%-6s comm_wait rank %d, box %d, %s.%s\n",
                "CREATE", sim->rank, i,
                PB_TYPENAME(type), PB_REQNAME(req));
        box_waitmpi(sim, &sim->box[i], type, dir, req);
    }
}
