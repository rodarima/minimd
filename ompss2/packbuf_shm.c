#define ENABLE_DEBUG 0
#include "types.h"
#include "log.h"
#include "safe.h"
#include "packbuf.h"

#include <string.h>
#include <stdlib.h>

void
packbuf_shm_init(PackBuf *pb, PackBuf *remote)
{
    pb->transport = PB_SHM;
    pb->shm.remote = remote;

    /* Allocate data with empty buf */
    pb->data = safe_calloc(1, sizeof(PackBufData));
    pb->nalloc = 0;

    packbuf_header_destroy_unsafe(pb);

    char tmp[256];
    strcpy(tmp, pb->name);

    /* Append custom shm info to the name */
    sprintf(pb->name, "%s[SHM remote=%p]", tmp, remote);
}

void
packbuf_shm_send(PackBuf *pb, enum pb_req reqtype)
{
    /* No need to do anything, as the header is already reset in
     * packbuf_send() */
    packbuf_switch(pb, PB_READY, PB_SENDING);
    packbuf_header_check(pb);
    packbuf_switch(pb, PB_SENDING, PB_READY);
}

static void
recv_natoms(PackBuf *pb)
{
    PackBuf *dst = pb;
    PackBuf *src = pb->shm.remote;

    size_t nbytes = packbuf_data_size(0, src->atomsize);
    memcpy(dst->data, src->data, nbytes);

    dst->natoms = src->natoms;
}

static void
recv_buf(PackBuf *pb)
{
    PackBuf *dst = pb;
    PackBuf *src = pb->shm.remote;

    dst->natoms = src->natoms;

    if (src->natoms != 0) {
        packbuf_grow(dst, src->natoms);
    }

    /* Always perform the copy, even with no atoms to transfer the header */
    size_t nbytes = packbuf_data_size(src->natoms, src->atomsize);
    memcpy(dst->data, src->data, nbytes);

    /* Set the natoms anyway */
    dst->natoms = src->natoms;
}

void
packbuf_shm_recv(PackBuf *pb, enum pb_req reqtype)
{
    PackBuf *dst = pb;
    PackBuf *src = pb->shm.remote;

    dbg("packbuf_shm_recv: natoms=%d reqtype=%d\n",
            src->natoms, reqtype);

    packbuf_switch(src, PB_READY, PB_COPYING);
    packbuf_switch(dst, PB_READY, PB_COPYING);

    if (reqtype == PB_NATOMS)
        recv_natoms(pb);
    else
        recv_buf(pb);

    /* We can already check the header here, as the transaction has just
     * finished */

    packbuf_header_check(pb);

    packbuf_switch(dst, PB_COPYING, PB_READY);
    packbuf_switch(src, PB_COPYING, PB_READY);
}
