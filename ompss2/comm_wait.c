#define ENABLE_DEBUG 0
#include "comm.h"
#include "types.h"
#include "log.h"
#include "packbuf.h"
#include "trace.h"

static void
box_waitmpi_task(Sim *sim, Box *box,
        enum pb_type type, enum pb_dir dir, enum pb_req req)
{
    dbg("%-6s comm_wait rank %d, box %d, %s.%s\n",
            "RUN", sim->rank, box->i,
            PB_TYPENAME(type), PB_REQNAME(req));

    packbuf_mpi_waitn(box->pb[type][dir], NNEIGH, req);

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
    for (int i = 0; i < NNEIGH; i++) {
        PackBuf *pb = box->pb[type][dir][i];
        #pragma oss task label("box_waitmpi set flag") \
            inout(pb->natoms) inout(pb->data) firstprivate(pb)
        {
            packbuf_debug_switch(pb, PB_READY, PB_WAITING);

            if (pb->trace[req].started) {
                pb->trace[req].t1 = MPI_Wtime();
                //trace_record(box->i * sim->nboxes + i, pb->trace[req].h,
                //        pb->trace[req].t0, pb->trace[req].t1, pb->trace[req].label, pb->trace[req].color);
                pb->trace[req].started = 0;
            }
            char label[1024];
            PackBufHeader *h = &pb->data->header;
            sprintf(label, "EVENT BEFORE box_waitmpi natoms=%d xnatoms=%d %s.%s.%s %s\n"
                    "header: magic=%d iter=%d srcbox=%d senddir=%d dstbox=%d icomm=%d\n"
		    "data: buf[0]=%e buf[1]=%e buf[2]=%e",
                pb->natoms, pb->data->xnatoms, PB_TYPENAME(type), PB_DIRNAME(dir), PB_REQNAME(req),
                pb->name,
                h->magic, h->iter, h->srcbox, h->senddir, h->dstbox, h->icomm,
		pb->data->buf[0], pb->data->buf[1], pb->data->buf[2]);
            trace_event(box->i * sim->nboxes + i, label, "#888800");
        }
    }

    #pragma oss task label("box_waitmpi") \
        inout({box->pb[type][dir][i]->natoms, i=0;NNEIGH}) \
        inout({box->pb[type][dir][i]->data,   i=0;NNEIGH})
    {
        //MPI_Barrier(MPI_COMM_WORLD);
        box_waitmpi_task(sim, box, type, dir, req);
    }

    for (int i = 0; i < NNEIGH; i++) {
        PackBuf *pb = box->pb[type][dir][i];
        #pragma oss task label("box_waitmpi set flag") \
            inout(pb->natoms) inout(pb->data) firstprivate(pb)
        {
            char label[1024];
            PackBufHeader *h = &pb->data->header;
            sprintf(label, "EVENT AFTER box_waitmpi natoms=%d xnatoms=%d %s.%s.%s %s\n"
                    "header: magic=%d iter=%d srcbox=%d senddir=%d dstbox=%d icomm=%d",
                pb->natoms, pb->data->xnatoms, PB_TYPENAME(type), PB_DIRNAME(dir), PB_REQNAME(req),
                pb->name,
                h->magic, h->iter, h->srcbox, h->senddir, h->dstbox, h->icomm);
            trace_event(box->i * sim->nboxes + i, label, "#880000");
        }
    }
}

void
comm_wait(Sim *sim,
        enum pb_type type, enum pb_dir dir, enum pb_req req)
{
    //dbg("rank%d -- comm_tidy -- recv recv_rvt buf\n", sim->rank);

    for (int i = 0; i < sim->nboxes; i++) {
        dbg("%-6s comm_wait rank %d, box %d, %s.%s.%s\n",
                "CREATE", sim->rank, i,
                PB_TYPENAME(type), PB_REQNAME(req), PB_DIRNAME(dir));
        box_waitmpi(sim, &sim->box[i], type, dir, req);
    }
}
