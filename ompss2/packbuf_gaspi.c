#define ENABLE_DEBUG 1
#include "log.h"
#include "types.h"
#include "packbuf.h"
#include "gaspi_check.h"

#include <GASPI.h>
#include <TAGASPI.h>

void
packbuf_gaspi_init(PackBuf *pb,
        double *newbuf,
        int sendseg, size_t send_offset_bytes,
        int recvseg, size_t recv_offset_bytes,
        size_t nalloc, int queue)
{
    pb->sendseg = sendseg;
    pb->recvseg = recvseg;
    pb->sendoffset = send_offset_bytes;
    pb->recvoffset = recv_offset_bytes;
    pb->buf = newbuf;
    pb->queue = queue;
    pb->gaspi = 1;
    pb->nalloc = nalloc;
}

static void
send_buf(PackBuf *pb)
{
    packbuf_switch(pb, PB_READY, PB_SENDING);

    dbg("packbuf_gaspi:send_buf: natoms=%d remoterank=%d tag=%d name='%s'\n",
            pb->natoms, pb->remoterank, pb->tag[PB_BUF], pb->name);

    if (pb->waitreq[PB_BUF])
        die("packbuf_gaspi_send_buf: buffer in use\n");

    /* Repeat until success */
    while (1) {
        gaspi_return_t ret = tagaspi_write_notify(
                pb->sendseg, pb->sendoffset,
                pb->remoterank,
                pb->recvseg, pb->recvoffset,
                pb->natoms * pb->atomsize,
                pb->tag[PB_BUF], 1,
                pb->queue);

        if (ret == GASPI_SUCCESS)
            break;

        if (ret != GASPI_QUEUE_FULL) {
            check_gaspi(ret, "tagaspi_write_notify",
                    __FILE__, __LINE__);
        }
    }

    /* FIXME: we may need to wait before overwriting the buffer */
    //pb->waitreq[PB_BUF] = 1;
    packbuf_switch(pb, PB_SENDING, PB_READY);
}

void
packbuf_gaspi_send(PackBuf *pb, enum pb_reqtype reqtype)
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
    packbuf_switch(pb, PB_READY, PB_RECVING);

    dbg("packbuf_gaspi_recv_buf: recvnatoms=%d remoterank=%d tag=%d name='%s'\n",
            pb->recvnatoms, pb->remoterank, pb->tag[PB_BUF], pb->name);

    if (pb->waitreq[PB_BUF])
        die("packbuf_gaspi_recv_buf: buffer in use\n");

    if (pb->recvnatoms < 0)
        die("bad recvnatoms\n");

    if (pb->recvnatoms > 0) {
        if (pb->recvnatoms > pb->nalloc)
            die("packbuf_gaspi_recv_buf: buffer of %d too small for %d atoms\n",
                    pb->nalloc, pb->recvnatoms);

        while (1) {
            gaspi_return_t ret = tagaspi_notify_async_wait(
                    pb->recvseg,
                    pb->tag[PB_BUF],
                    GASPI_NOTIFICATION_IGNORE);

            if (ret == GASPI_SUCCESS)
                break;

            if (ret != GASPI_QUEUE_FULL) {
                check_gaspi(ret, "tagaspi_notify_async_wait",
                        __FILE__, __LINE__);
            }
        }

        //pb->waitreq[PB_BUF] = 1;
    }

    pb->natoms = pb->recvnatoms;
    packbuf_switch(pb, PB_RECVING, PB_READY);
}

void
packbuf_gaspi_recv(PackBuf *pb, enum pb_reqtype reqtype)
{
    if (reqtype == PB_NATOMS) {
        die("not implemented\n");
    } else {
        recv_buf(pb);
    }
}
