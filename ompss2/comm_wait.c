#define ENABLE_DEBUG 0
#include "types.h"
#include "log.h"
#include "packbuf.h"

static void
box_waitmpi_task(Sim *sim, Box *box, enum pb_type type, enum pb_reqtype reqtype)
{
    dbg("%-6s comm_wait rank %d, box %d, %s.%s buf'=%p\n",
            "RUN", sim->rank, box->i,
            PB_TYPENAME(type), PB_REQTYPENAME(reqtype),
            &box->pb[type][0]->buf);

    for (int i = 0; i < NNEIGH; i++) {
        PackBuf *pb = box->pb[type][i];
        packbuf_debug_switch(pb, PB_READY, PB_WAITING);
    }

    packbuf_mpi_waitn(box->pb[type], NNEIGH, reqtype);

    dbg("%-6s comm_wait rank %d, box %d, %s.%s\n",
            "DONE", sim->rank, box->i,
            PB_TYPENAME(type), PB_REQTYPENAME(reqtype));


    for (int i = 0; i < NNEIGH; i++) {
        PackBuf *pb = box->pb[type][i];
        packbuf_debug_switch(pb, PB_WAITING, PB_READY);
    }
}

static void
box_waitmpi(Sim *sim, Box *box, enum pb_type type, enum pb_reqtype reqtype)
{
    if (!NEED_EXPLICIT_WAIT)
        return;

    #pragma oss task label("box_waitmpi") \
        inout({box->pb[type][i]->natoms, i=0;NNEIGH}) \
        inout({box->pb[type][i]->buf,    i=0;NNEIGH})
    {
        box_waitmpi_task(sim, box, type, reqtype);
    }
}

void
comm_wait(Sim *sim, enum pb_type type, enum pb_reqtype reqtype)
{
    //dbg("rank%d -- comm_tidy -- recv recv_rvt buf\n", sim->rank);

    for (int i = 0; i < sim->nboxes; i++) {
        dbg("%-6s comm_wait rank %d, box %d, %s.%s buf=%p\n",
                "CREATE", sim->rank, i,
                PB_TYPENAME(type), PB_REQTYPENAME(reqtype),
                &sim->box[i].pb[type][0]->buf);
        box_waitmpi(sim, &sim->box[i], type, reqtype);
    }
}
