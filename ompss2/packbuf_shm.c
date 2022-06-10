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

    char tmp[256];
    strcpy(tmp, pb->name);

    /* Append custom shm info to the name */
    sprintf(pb->name, "%s[SHM remote=%p]", tmp, remote);
}

void
packbuf_shm_send(PackBuf *pb, enum pb_req reqtype)
{
    die("%s: use only recv with SHM\n", pb->name);
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
        die("%s: non-sense receive only natoms via SHM\n", pb->name);

    if (src->natoms != 0) {
        packbuf_grow(dst, src->natoms);
    }

    /* Always perform the copy, even with no atoms to transfer the header */
    size_t nbytes = packbuf_data_size(src->natoms, src->atomsize);
    memcpy(dst->data, src->data, nbytes);

    dst->natoms = src->natoms;

    packbuf_switch(dst, PB_COPYING, PB_READY);
    packbuf_switch(src, PB_COPYING, PB_READY);
}
