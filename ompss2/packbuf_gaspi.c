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

void
packbuf_gaspi_send_buf(PackBuf *pb)
{
    packbuf_switch(pb, PB_READY, PB_SENDING);

    dbg("packbuf_gaspi_send_buf: natoms=%d remoterank=%d tag=%d\n",
            pb->natoms, pb->remoterank, pb->tag);

    if (pb->waitreq)
        die("packbuf_gaspi_send_buf: buffer in use\n");

    /* Repeat until success */
    while (1) {
        gaspi_return_t ret = tagaspi_write_notify(
                pb->sendseg, pb->sendoffset,
                pb->remoterank,
                pb->recvseg, pb->recvoffset,
                pb->natoms * pb->atomsize,
                pb->tag, 1,
                pb->queue);

        if (ret == GASPI_SUCCESS)
            break;

        if (ret != GASPI_QUEUE_FULL) {
            check_gaspi(ret, "tagaspi_write_notify",
                    __FILE__, __LINE__);
        }
    }

    /* Always set the waitreq flag with tagaspi */
    pb->waitreq = 1;
    packbuf_switch(pb, PB_SENDING, PB_READY);
}

void
packbuf_gaspi_recv_buf(PackBuf *pb, int natoms)
{
    packbuf_switch(pb, PB_READY, PB_RECVING);

    dbg("packbuf_gaspi_recv_buf: natoms=%d remoterank=%d tag=%d\n",
            natoms, pb->remoterank, pb->tag);

    if (pb->waitreq)
        die("packbuf_gaspi_recv_buf: buffer in use\n");

    if (natoms > 0) {
        if (natoms < pb->nalloc)
            die("packbuf_gaspi_recv_buf: buffer too small for %d atoms\n",
                    natoms);

        while (1) {
            gaspi_return_t ret = tagaspi_notify_async_wait(
                    pb->recvseg,
                    pb->tag,
                    GASPI_NOTIFICATION_IGNORE);

            if (ret == GASPI_SUCCESS)
                break;

            if (ret != GASPI_QUEUE_FULL) {
                check_gaspi(ret, "tagaspi_notify_async_wait",
                        __FILE__, __LINE__);
            }
        }

        pb->waitreq = 1;
    }

    pb->natoms = natoms;
    packbuf_switch(pb, PB_RECVING, PB_READY);
}
