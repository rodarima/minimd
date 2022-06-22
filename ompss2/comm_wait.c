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

    packbuf_waitn(box->pb[type][dir], NNEIGH, req);

    for (int i = 0; i < NNEIGH; i++) {
        PackBuf *pb = box->pb[type][dir][i];
        packbuf_debug_switch(pb, PB_WAITING, PB_READY);
    }

    dbg("%-6s comm_wait rank %d, box %d, %s.%s\n",
            "DONE", sim->rank, box->i,
            PB_TYPENAME(type), PB_REQNAME(req));
}

static void
box_waitmpi(Sim *sim, Box *box,
        enum pb_type type, enum pb_dir dir, enum pb_req req)
{
    /* Set the waiting state as soon as we can in each PackBuf
     * independently, so we can look the state in the debugger */
    for (int i = 0; i < NNEIGH; i++) {
        PackBuf *pb = box->pb[type][dir][i];
        #pragma oss task label("box_waitmpi set flag") \
            inout(pb->natoms) inout(pb->data) firstprivate(pb)
        {
            packbuf_debug_switch(pb, PB_READY, PB_WAITING);
        }
    }

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
    for (int i = 0; i < sim->nboxes; i++)
        box_waitmpi(sim, &sim->box[i], type, dir, req);
}
