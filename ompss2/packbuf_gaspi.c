#define ENABLE_DEBUG 1
#include "log.h"
#include "types.h"
#include "packbuf.h"
#include "gaspi_check.h"
#include "trace.h"

#include <GASPI.h>
#include <TAGASPI.h>
#include <string.h>

static int
check_tag(int tag)
{
    gaspi_number_t maxtag;
    CHECK(gaspi_notification_num(&maxtag));

    if (tag >= (int) maxtag)
        die("GASPI tag exceed limit: %d >= %d\n", tag, maxtag);

    return tag;
}

void
packbuf_gaspi_init(PackBuf *pb,
        PackBufData *newdata,
        int sendseg, size_t send_offset_bytes,
        int recvseg, size_t recv_offset_bytes,
        size_t nalloc, int queue, int tag)
{
    PackBufGASPI *pbg = &pb->gaspi;
    pbg->sendseg = sendseg;
    pbg->recvseg = recvseg;
    pbg->sendoffset = send_offset_bytes;
    pbg->recvoffset = recv_offset_bytes;
    pbg->queue = queue;
    pbg->nid = check_tag(tag);

    pb->nalloc = nalloc;
    pb->data = newdata;
    pb->transport = PB_GASPI;

    /* Append custom shm info to the name */
    char tmp[1024];
    strcpy(tmp, pb->name);
    sprintf(pb->name, "%s[GASPI sendseg=%d recvseg=%d "
            "sendoffset=%lu recvoffset=%lu "
            "queue=%d nid=%d]", tmp,
            pbg->sendseg, pbg->recvseg,
            pbg->sendoffset, pbg->recvoffset,
            pbg->queue, pbg->nid);
}

static void
send_buf(PackBuf *pb)
{
    PackBufGASPI *pbg = &pb->gaspi;
    packbuf_switch(pb, PB_READY, PB_SENDING);

//    dbg("packbuf_gaspi:send_buf: natoms=%d remoterank=%d nid=%d sendoffset=%lu recvoffset=%lu name='%s'\n",
//            pb->natoms, pb->remoterank, pbg->nid,
//            pbg->sendoffset, pbg->recvoffset, pb->name);
    int bytes = sizeof(*pb->data) +
        pb->natoms * pb->atomsize * sizeof(double);

    /* Repeat until success */
    while (1) {
        PackBufHeader *h = &pb->data->header;
        char label[1024];
        sprintf(label, "tagaspi_write_notify(sendseg=%d, sendoffset=%lu, remoterank=%d, \n"
                "recvseg=%d, recvoffset=%lu, nbytes=%d, nid=%d, one=%d, queue=%d) ENTER\n"
                "header: magic=%d iter=%d srcbox=%d senddir=%d dstbox=%d icomm=%d\n"
                "magic=%d natoms=%d buf[0]=%e buf[1]=%e buf[1]=%e",
                pbg->sendseg, pbg->sendoffset,
                pb->remoterank,
                pbg->recvseg, pbg->recvoffset,
                bytes,
                pbg->nid, 1,
                pbg->queue,
                h->magic, h->iter, h->srcbox, h->senddir, h->dstbox, h->icomm,
                pb->data->header.magic,
                pb->natoms,
                pb->data->buf[0], pb->data->buf[1], pb->data->buf[2]);

        trace_event(pb->box * NNEIGH + pb->ineigh, label, "#33ff00");

		dbg("tagaspi_write_notify(localseg=%d, localoffset=%lu, remoterank=%d, \n"
				"  remoteseg=%d, remoteoffset=%lu, nbytes=%d, nid=%d, one=%d, queue=%d) ENTER\n"
				"  header: magic=%d iter=%d srcbox=%d senddir=%d dstbox=%d icomm=%d\n"
				"  magic=%d natoms=%d buf[0]=%e buf[1]=%e buf[1]=%e\n",
				pbg->sendseg, pbg->sendoffset,
				pb->remoterank,
				pbg->recvseg, pbg->recvoffset,
				bytes,
				pbg->nid, 1,
				pbg->queue,
				h->magic, h->iter, h->srcbox, h->senddir, h->dstbox, h->icomm,
				pb->data->header.magic,
				pb->natoms,
				pb->data->buf[0], pb->data->buf[1], pb->data->buf[2]);

        gaspi_return_t ret = tagaspi_write_notify(
                pbg->sendseg, pbg->sendoffset,
                pb->remoterank,
                pbg->recvseg, pbg->recvoffset,
                bytes,
                pbg->nid, 1,
                pbg->queue);

        sprintf(label, "tagaspi_write_notify(sendseg=%d, sendoffset=%lu, remoterank=%d, \n"
                "recvseg=%d, recvoffset=%lu, count=%d, nid=%d, one=%d, queue=%d) EXIT=%d",
                pbg->sendseg, pbg->sendoffset,
                pb->remoterank,
                pbg->recvseg, pbg->recvoffset,
                bytes,
                pbg->nid, 1,
                pbg->queue, ret);

        trace_event(pb->box * NNEIGH + pb->ineigh, label, "#00ff00");

        if (ret == GASPI_SUCCESS)
            break;

        if (ret != GASPI_QUEUE_FULL) {
            check_gaspi(ret, "tagaspi_write_notify",
                    __FILE__, __LINE__);
        }
    }

    /* TODO: Set in_transfer and wait before attempting another send
     * until the remote buffer has been emptied. */

    packbuf_switch(pb, PB_SENDING, PB_READY);
}

void
packbuf_gaspi_send(PackBuf *pb, enum pb_req reqtype)
{
    if (reqtype == PB_NATOMS) {
        die("not implemented\n");
    } else {
        send_buf(pb);
    }
}

static void
recv_buf(PackBuf *pb)
{
    PackBufGASPI *pbg = &pb->gaspi;
    packbuf_switch(pb, PB_READY, PB_RECVING);

    /* xnatoms will be overwritten by the message */
    int recvnatoms = pb->natoms;

    dbg("packbuf_gaspi_recv_buf: recvnatoms=%d remoterank=%d nid=%d sendoffset=%lu recvoffset=%lu name='%s'\n",
            pb->natoms, pb->remoterank, pbg->nid,
            pbg->sendoffset, pbg->recvoffset, pb->name);

    if (recvnatoms < 0)
        die("%s: negative recvnatoms=%d\n", pb->name, recvnatoms);

    if (recvnatoms > 0) {
        if (recvnatoms > pb->nalloc) {
            die("packbuf_gaspi_recv_buf: buffer of %d too small for %d atoms\n",
                    pb->nalloc, recvnatoms);
        }

        while (1) {

            char label[1024];
            sprintf(label, "tagaspi_notify_async_wait(seg=%d, nid=%d) ENTER\n"
                    "magic=%d buf[0]=%e buf[1]=%e buf[2]=%e",
                    pbg->recvseg, pbg->nid,
                    pb->data->header.magic,
                    pb->data->buf[0], pb->data->buf[1], pb->data->buf[2]);
            trace_event(pb->box * NNEIGH + pb->ineigh, label, "#00ffff");

            dbg("tagaspi_notify_async_wait(seg=%d, nid=%d)\n"
                    "  magic=%d buf[0]=%e buf[1]=%e buf[2]=%e\n",
                    pbg->recvseg, pbg->nid,
                    pb->data->header.magic,
                    pb->data->buf[0], pb->data->buf[1], pb->data->buf[2]);

            gaspi_return_t ret = tagaspi_notify_async_wait(
                    pbg->recvseg,
                    pbg->nid,
                    GASPI_NOTIFICATION_IGNORE);

            sprintf(label, "tagaspi_notify_async_wait(seg=%d, nid=%d) magic=%d EXIT=%d",
                    pbg->recvseg, pbg->nid, pb->data->header.magic, ret);
            trace_event(pb->box * NNEIGH + pb->ineigh, label, "#00ffff");

            if (ret == GASPI_SUCCESS)
                break;

            if (ret != GASPI_QUEUE_FULL) {
                check_gaspi(ret, "tagaspi_notify_async_wait",
                        __FILE__, __LINE__);
            }
        }

        pb->in_transfer[PB_BUF] = 1;
    }

    /* We cannot check the header yet as only after the task is released we
     * would have received the data */

    packbuf_switch(pb, PB_RECVING, PB_READY);
}

void
packbuf_gaspi_recv(PackBuf *pb, enum pb_req reqtype)
{
    if (reqtype == PB_NATOMS) {
        die("not implemented\n");
    } else {
        recv_buf(pb);
    }
}

void
packbuf_gaspi_waitn(PackBuf **pbs, int n, enum pb_req req)
{
    for (int i = 0; i < n; i++) {
        PackBuf *pb = pbs[i];

        if (pb->transport != PB_GASPI || !pb->in_transfer[req])
            continue;

        packbuf_switch(pb, PB_READY, PB_WAITING);

        if (req != PB_BUF)
            die("packbuf_gaspi_waitn: only PB_BUF supported\n");

        /* Ensure the header is sane */
        if (pb->dir == PB_RECV) {
            // FIXME: disabled for testing
            //packbuf_check_header(pb, -666);

            /* No need to set natoms as we only exchange PB_BUF */
            //pb->natoms = pb->data->xnatoms;
        }

        pb->in_transfer[req] = 0;

        packbuf_switch(pb, PB_WAITING, PB_READY);
    }
}
