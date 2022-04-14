#define ENABLE_DEBUG 1
#include "types.h"
#include "log.h"
#include "packbuf.h"

static void
box_waitmpi(Sim *sim, Box *box, enum pb_type type, enum pb_reqtype reqtype)
{
    if (!NEED_EXPLICIT_WAIT)
        return;

    /* FIXME ! */
    #pragma oss taskwait

    if (type == PB_SEND_R || type == PB_SEND_RT || type == PB_SEND_RVT) {
        #pragma oss task label("box_waitmpi_send") \
            in({box->pb[type][i]->natoms, i=0;NNEIGH}) \
            in({box->pb[type][i]->buf,    i=0;NNEIGH})
        {
            packbuf_mpi_waitn(box->pb[type], NNEIGH, reqtype);
        }
    } else {
        #pragma oss task label("box_waitmpi_recv") \
            inout({box->pb[type][i]->natoms, i=0;NNEIGH}) \
            inout({box->pb[type][i]->buf,    i=0;NNEIGH})
        {
            packbuf_mpi_waitn(box->pb[type], NNEIGH, reqtype);
        }
    }

    /* FIXME ! */
    #pragma oss taskwait
}

void
comm_wait(Sim *sim, enum pb_type type, enum pb_reqtype reqtype)
{
    //dbg("rank%d -- comm_tidy -- recv recv_rvt buf\n", sim->rank);

    for (int i = 0; i < sim->nboxes; i++)
        box_waitmpi(sim, &sim->box[i], type, reqtype);
}
