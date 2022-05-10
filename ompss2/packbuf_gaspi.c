#define ENABLE_DEBUG 0
#include "log.h"
#include "types.h"
#include "packbuf.h"
#include "gaspi_check.h"

#include <GASPI.h>
#include <TAGASPI.h>

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
    pb->mode = PB_GASPI;
}

static void
send_buf(PackBuf *pb)
{
    PackBufGASPI *pbg = &pb->gaspi;
    packbuf_switch(pb, PB_READY, PB_SENDING);

    dbg("packbuf_gaspi:send_buf: natoms=%d remoterank=%d nid=%d name='%s'\n",
            pb->natoms, pb->remoterank, pbg->nid, pb->name);

    if (pb->waitreq[PB_BUF])
        die("packbuf_gaspi_send_buf: buffer in use\n");

    /* Repeat until success */
    while (1) {
        gaspi_return_t ret = tagaspi_write_notify(
                pbg->sendseg, pbg->sendoffset,
                pb->remoterank,
                pbg->recvseg, pbg->recvoffset,
                pb->natoms * pb->atomsize,
                pbg->nid, 1,
                pbg->queue);

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

    dbg("packbuf_gaspi_recv_buf: recvnatoms=%d remoterank=%d nid=%d name='%s'\n",
            recvnatoms, pb->remoterank, pbg->nid, pb->name);

    if (pb->waitreq[PB_BUF])
        die("packbuf_gaspi_recv_buf: buffer in use\n");

    if (recvnatoms < 0)
        die("%s: negative recvnatoms=%d\n", pb->name, recvnatoms);

    if (recvnatoms > 0) {
        if (recvnatoms > pb->nalloc)
            die("packbuf_gaspi_recv_buf: buffer of %d too small for %d atoms\n",
                    pb->nalloc, recvnatoms);

        while (1) {
            gaspi_return_t ret = tagaspi_notify_async_wait(
                    pbg->recvseg,
                    pbg->nid,
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
