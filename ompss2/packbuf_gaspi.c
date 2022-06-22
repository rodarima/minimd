#define ENABLE_DEBUG 0
#include "log.h"
#include "types.h"
#include "packbuf.h"
#include "gaspi_check.h"

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
        int local_seg, size_t local_offset_bytes,
        int remote_seg, size_t remote_offset_bytes,
        size_t nalloc, int queue, int tag)
{
    PackBufGASPI *pbg = &pb->gaspi;
    pbg->local_seg = local_seg;
    pbg->local_offset = local_offset_bytes;
    pbg->remote_seg = remote_seg;
    pbg->remote_offset = remote_offset_bytes;
    pbg->queue = queue;
    pbg->nid = check_tag(tag);

    if (pb->dir == PB_SEND) {
        pbg->src_seg = local_seg;
        pbg->src_offset = local_offset_bytes;
        pbg->dst_seg = remote_seg;
        pbg->dst_offset = remote_offset_bytes;
    } else {
        pbg->src_seg = remote_seg;
        pbg->src_offset = remote_offset_bytes;
        pbg->dst_seg = local_seg;
        pbg->dst_offset = local_offset_bytes;
    }

    pb->nalloc = nalloc;
    pb->data = newdata;
    pb->transport = PB_GASPI;

    packbuf_header_destroy_unsafe(pb);

    /* Append custom shm info to the name */
    char tmp[1024];
    strcpy(tmp, pb->name);
    sprintf(pb->name, "%s[GASPI local_seg=%d remote_seg=%d "
            "local_offset=%lu remote_offset=%lu "
            "queue=%d nid=%d]", tmp,
            pbg->local_seg, pbg->remote_seg,
            pbg->local_offset, pbg->remote_offset,
            pbg->queue, pbg->nid);
}

static void
send_buf(PackBuf *pb)
{
    PackBufGASPI *pbg = &pb->gaspi;
    packbuf_switch(pb, PB_READY, PB_SENDING);

    int bytes = sizeof(*pb->data) +
        pb->natoms * pb->atomsize * sizeof(double);

    /* Repeat until success */
    while (1) {
        PackBufHeader *h = &pb->data->header;

        gaspi_return_t ret = tagaspi_write_notify(
                pbg->local_seg, pbg->local_offset,
                pb->remote.rank,
                pbg->remote_seg, pbg->remote_offset,
                bytes,
                pbg->nid, 1,
                pbg->queue);

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

    dbg("packbuf_gaspi_recv_buf: recvnatoms=%d name=%s\n",
            recvnatoms, pb->name);

    if (recvnatoms < 0)
        die("%s: negative recvnatoms=%d\n", pb->name, recvnatoms);

    if (recvnatoms > 0) {
        if (recvnatoms > pb->nalloc) {
            die("packbuf_gaspi_recv_buf: buffer of %d too small for %d atoms\n",
                    pb->nalloc, recvnatoms);
        }

        while (1) {

            dbg("tagaspi_notify_async_wait(seg=%d, nid=%d)\n"
                    "  magic=%d buf[0]=%e buf[1]=%e buf[2]=%e\n",
                    pbg->local_seg, pbg->nid,
                    pb->data->header.magic,
                    pb->data->buf[0], pb->data->buf[1], pb->data->buf[2]);

            gaspi_return_t ret = tagaspi_notify_async_wait(
                    pbg->local_seg,
                    pbg->nid,
                    GASPI_NOTIFICATION_IGNORE);

            if (ret == GASPI_SUCCESS)
                break;

            if (ret != GASPI_QUEUE_FULL) {
                check_gaspi(ret, "tagaspi_notify_async_wait",
                        __FILE__, __LINE__);
            }
        }

        pb->in_transfer[PB_BUF] = 1;
    }

    /* We cannot check the header yet as only after the task is released
     * we would have received the data */

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
packbuf_gaspi_signal(PackBuf *pb)
{
    PackBufGASPI *pbg = &pb->gaspi;
    packbuf_switch(pb, PB_READY, PB_SENDING);

    /* Only valid for recv PackBuf */
    if (pb->dir != PB_RECV)
        die("packbuf_gaspi_signal: only valid for PB_RECV\n");

    /* Repeat until success */
    while (1) {
        /* This is counter-intuitive, because we are sending a message
         * to the source endpoint to acknowledge that the receiving end
         * is ready to receive more data. */
        gaspi_return_t ret = tagaspi_notify(
                pbg->src_seg, pb->remote.rank,
                pbg->nid, 1, pbg->queue);

        if (ret == GASPI_SUCCESS)
            break;

        if (ret != GASPI_QUEUE_FULL) {
            check_gaspi(ret, "tagaspi_notify", __FILE__, __LINE__);
        }
    }

    packbuf_switch(pb, PB_SENDING, PB_READY);
}

void
packbuf_gaspi_linger(PackBuf *pb)
{
    PackBufGASPI *pbg = &pb->gaspi;
    packbuf_switch(pb, PB_READY, PB_WAITING);

    if (pb->dir != PB_SEND)
        die("packbuf_gaspi_linger: only valid for PB_SEND\n");

    while (1) {

        /* FIXME: Ensure we don't receive crossed notifications from
         * other parts of the simulation */

        /* Wait for the notification in the source segment */
        gaspi_return_t ret = tagaspi_notify_async_wait(
                pbg->src_seg, pbg->nid,
                GASPI_NOTIFICATION_IGNORE);

        if (ret == GASPI_SUCCESS)
            break;

        if (ret != GASPI_QUEUE_FULL) {
            check_gaspi(ret, "tagaspi_notify_async_wait",
                    __FILE__, __LINE__);
        }
    }

    packbuf_switch(pb, PB_WAITING, PB_READY);
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
            packbuf_header_check(pb);

            /* No need to set natoms as we only exchange PB_BUF */
            //pb->natoms = pb->data->xnatoms;
        }

        pb->in_transfer[req] = 0;

        packbuf_switch(pb, PB_WAITING, PB_READY);
    }
}
