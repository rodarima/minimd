#define ENABLE_DEBUG 0
#include "types.h"
#include "log.h"
#include "safe.h"
#include "packbuf.h"

#include <string.h>
#include <stdlib.h>
#include <stddef.h>

#define PACKBUF_INCR 2000

void
packbuf_switch(PackBuf *pb, enum packbuf_state prev, enum packbuf_state next)
{
    if (ENABLE_PACKBUF_STATE) {
        if (pb->state != prev) {
            die("packbuf in state %d, expected %d (switching to %d)\n",
                    pb->state, prev, next);
        }
        pb->state = next;
    }
}

void
packbuf_debug_switch(PackBuf *pb, enum packbuf_state prev, enum packbuf_state next)
{
    if (ENABLE_PACKBUF_DEBUG_STATE) {
        if (pb->debug_state != prev) {
            die("packbuf in debug_state %d, expected %d (switching to %d)\n",
                    pb->debug_state, prev, next);
        }

        pb->debug_state = next;
    }
}

size_t
packbuf_data_size(size_t natoms, size_t atomdoubles)
{
    size_t datasize = natoms * atomdoubles * sizeof(double);
    return sizeof(PackBufData) + datasize;
}

/* Grows the buffer so that the allocated capacity can hold at least n
 * atoms. */
void
packbuf_grow(PackBuf *pb, int n)
{
    if (n < pb->natoms)
        n = pb->natoms;

    if (pb->nalloc < n) {
        if (pb->transport == PB_GASPI)
            die("packbuf_grow cannot grow GASPI buffer\n");

        if (n == 0)
            die("cannot allocate 0 atoms\n");

        //if (pb->nalloc + PACKBUF_INCR < n)
        //    n = pb->nalloc + PACKBUF_INCR;

        size_t newalloc = packbuf_data_size(n, pb->atomsize);

        pb->data = safe_realloc(pb->data, newalloc);

        /* Also grow the selection buffer if enabled */
        if (pb->enable_sel) {
            pb->sel = (int *) safe_realloc(pb->sel, sizeof(int) * n);
        }

        pb->nalloc = n;
    }
}

/* Ensures the buffer can hold at least nextra additional atoms */
static void
packbuf_grow_extra(PackBuf *pb, int nextra)
{
    packbuf_grow(pb, pb->natoms + nextra);
}

/* FIXME: move to .h so the compiler can optimize the constant NULL
 * pointers */
void
packbuf_add(PackBuf *pb, Vec *r, Vec *v, int *type)
{
    packbuf_switch(pb, PB_READY, PB_ADDING);
    /* Ensure we have room for another atom */
    packbuf_grow_extra(pb, 1);

    if (pb->in_transfer[PB_BUF] || pb->in_transfer[PB_NATOMS])
        die("packbuf_add: buffer in use\n");

    int j = pb->natoms * pb->atomsize;

    if (r != NULL) {
        for (int d = X; d <= Z; d++)
            pb->data->buf[j++] = r[0][d];
    }

    if (v != NULL) {
        for (int d = X; d <= Z; d++)
            pb->data->buf[j++] = v[0][d];
    }

    if (type != NULL) {
        /* FIXME: We are sending the type as a double */
        pb->data->buf[j++] = (double) type[0];
    }

    pb->natoms++;

//    dbg("packbuf %p now has %d atoms (alloc %d)\n",
//            pb, pb->natoms, pb->nalloc);

    if (j != pb->natoms * pb->atomsize)
        die("packbuf_add atom size mismatch\n");

    packbuf_switch(pb, PB_ADDING, PB_READY);
}

void
packbuf_add_sel(PackBuf *pb, Vec *r, Vec *v, int *type, int iatom)
{
    int ncur = pb->natoms;
    packbuf_add(pb, r, v, type);
    pb->sel[ncur] = iatom;
}

void
packbuf_unpack(PackBuf *pb, Vec *r, Vec *v, int *types)
{
    packbuf_switch(pb, PB_READY, PB_UNPACKING);

    if (pb->in_transfer[PB_BUF] || pb->in_transfer[PB_NATOMS])
        die("packbuf_unpack: buffer in use, in_transfer(buf=%d natoms=%d) %s\n",
                pb->in_transfer[PB_BUF], pb->in_transfer[PB_NATOMS],
                pb->name);

    for (int i = 0, j = 0; i < pb->natoms; i++) {
        if (r != NULL) {
            for (int d = X; d <= Z; d++)
                r[i][d] = pb->data->buf[j++];
        }

        if (v != NULL) {
            for (int d = X; d <= Z; d++)
                v[i][d] = pb->data->buf[j++];
        }

        if (types != NULL) {
            types[i] = (int) pb->data->buf[j++];
        }
    }
    packbuf_switch(pb, PB_UNPACKING, PB_READY);
}

void
packbuf_unpack_sel(PackBuf *pb, Vec *r, Vec *v, int *types, int *sel)
{
    packbuf_switch(pb, PB_READY, PB_UNPACKING);

    if (pb->in_transfer[PB_BUF] || pb->in_transfer[PB_NATOMS])
        die("packbuf_unpack_sel: buffer in use\n");

    for (int i = 0, j = 0; i < pb->natoms; i++) {
        if (r != NULL) {
            for (int d = X; d <= Z; d++)
                r[sel[i]][d] = pb->data->buf[j++];
        }

        if (v != NULL) {
            for (int d = X; d <= Z; d++)
                v[sel[i]][d] = pb->data->buf[j++];
        }

        if (types != NULL) {
            types[sel[i]] = (int) pb->data->buf[j++];
        }
    }
    packbuf_switch(pb, PB_UNPACKING, PB_READY);
}

void
packbuf_clear(PackBuf *pb)
{
    packbuf_switch(pb, PB_READY, PB_CLEANING);

    if (pb->in_transfer[PB_BUF] || pb->in_transfer[PB_NATOMS])
        die("packbuf_unpack_sel: buffer in use\n");

    pb->natoms = 0;

    packbuf_switch(pb, PB_CLEANING, PB_READY);
}

void
packbuf_header_reset(PackBuf *pb)
{
    PackBufHeader *h = &pb->data->header;

    h->magic = PB_MAGIC_OK;
    h->iseq = pb->iseq;

    memcpy(&h->src, pb->src, sizeof(Endpoint));
    memcpy(&h->dst, pb->dst, sizeof(Endpoint));

    /* We may want to include transport specific data in the header as
     * well */
}

void
packbuf_header_destroy_unsafe(PackBuf *pb)
{
    PackBufHeader *h = &pb->data->header;

    /* Destroy the whole thing */
    memset(h, 0xff, sizeof(*h));

    h->magic = PB_MAGIC_DESTROYED;
}

void
packbuf_header_destroy(PackBuf *pb)
{
    PackBufHeader *h = &pb->data->header;

    if (pb->in_transfer[PB_BUF] || pb->in_transfer[PB_NATOMS])
        die("cannot destroy header: in transfer is set %s", pb->name);

    if (h->magic != PB_MAGIC_OK) {
        die("cannot destroy header: magic is not OK (%d) %s\n",
                h->magic, pb->name);
    }

    packbuf_header_destroy_unsafe(pb);
}

/** Ensure the header matches with the expected values */
void
packbuf_header_check(PackBuf *pb)
{
    PackBufHeader *h = &pb->data->header;

    if (h->magic != PB_MAGIC_OK) {
        die("%s wrong magic %d (expected %d)\n",
                pb->name, h->magic, PB_MAGIC_OK);
    }

    if (!endpoint_is_same(&h->src, pb->src)) {
        die("%s source endpoint mismatch:\n"
                "  header:   %s\n"
                "  expected: %s\n",
                pb->name, h->src.name, pb->src->name);
    }

    if (!endpoint_is_same(&h->dst, pb->dst)) {
        die("%s destination endpoint mismatch:\n"
                "  header:   %s\n"
                "  expected: %s\n",
                pb->name, h->dst.name, pb->dst->name);
    }

    if (h->iseq != pb->iseq) {
        die("%s iseq number mismatch: header %d, expected %d\n",
                pb->name, h->iseq, pb->iseq);
    }

    dbg("%s header ok\n", pb->name);
}

void
packbuf_send(PackBuf *pb, enum pb_req reqtype)
{
    if (pb->in_transfer[reqtype])
        die("packbuf_send: already transferring data %s\n", pb->name);

    pb->iseq++;
    dbg("packbuf_send: increased iseq to %d for %s\n",
            pb->iseq, pb->name);

    packbuf_header_reset(pb);

    switch (pb->transport) {
        case PB_GASPI: packbuf_gaspi_send(pb, reqtype); break;
        case PB_MPI: packbuf_mpi_send(pb, reqtype); break;
        case PB_SHM: packbuf_shm_send(pb, reqtype); break;
        default: die("packbuf_send: bad transport\n");
    }
}

void
packbuf_recv(PackBuf *pb, enum pb_req reqtype)
{
    if (pb->in_transfer[reqtype])
        die("packbuf_recv: already transferring data %s\n", pb->name);

    pb->iseq++;
    dbg("packbuf_recv: increased iseq to %d for %s\n",
            pb->iseq, pb->name);

    /* We cannot destroy the header here, as it may be overwritten
     * remotely by GASPI */

    switch (pb->transport) {
        case PB_GASPI: packbuf_gaspi_recv(pb, reqtype); break;
        case PB_MPI: packbuf_mpi_recv(pb, reqtype); break;
        case PB_SHM: packbuf_shm_recv(pb, reqtype); break;
        default: die("packbuf_recv: bad transport\n");
    }
}

void
packbuf_init(PackBuf *pb, int enable_sel, int atomsize,
        Endpoint *local, Endpoint *remote)
{
    memset(pb, 0, sizeof(*pb));

    pb->atomsize = atomsize;
    pb->enable_sel = enable_sel;
    pb->transport = PB_BAD;
    pb->dir = local->dir;
    pb->type = local->type;
    pb->data = NULL;

    memcpy(&pb->local, local, sizeof(Endpoint));
    memcpy(&pb->remote, remote, sizeof(Endpoint));

    /* Setup source/destination aliases */
    if (pb->dir == PB_SEND) {
        pb->src = &pb->local;
        pb->dst = &pb->remote;
    } else {
        pb->dst = &pb->local;
        pb->src = &pb->remote;
    }

    sprintf(pb->name, "PackBuf{type=%s dir=%s "
            "local={rank=%d box=%d index=%d type=%s dir=%s} "
            "remote={rank=%d box=%d index=%d type=%s dir=%s}}",
            PB_TYPENAME(pb->type), PB_DIRNAME(pb->dir),
            local->rank, local->box, local->index,
            PB_TYPENAME(local->type), PB_DIRNAME(local->dir),
            remote->rank, remote->box, remote->index,
            PB_TYPENAME(remote->type), PB_DIRNAME(remote->dir));

    packbuf_switch(pb, PB_GARBAGE, PB_READY);
}

void
packbuf_waitn(PackBuf **pbs, int n, enum pb_req reqtype)
{
    /* The array of pbs may contain mixed transports, so we just call
     * each implementation sequentially. Each transport will only
     * operate on its own transport PackBuf buffers */
    packbuf_mpi_waitn(pbs, n, reqtype);
    packbuf_gaspi_waitn(pbs, n, reqtype);

    /* Ensure all in_transfer flags are cleared */
    for (int i = 0; i < n; i++)
        if (pbs[i]->in_transfer[reqtype])
            die("packbuf_waitn: still in transfer %s\n", pbs[i]->name);
}

/** Wait until the remote buffer can be written */
void
packbuf_linger(PackBuf *pb)
{
    switch (pb->transport) {
        case PB_GASPI: packbuf_gaspi_linger(pb); break;
        default: break;
    }
}

/** Signals that the PackBuf is ready for receiving data. */
void
packbuf_signal(PackBuf *pb)
{
    if (pb->dir != PB_RECV)
        die("packbuf_signal: should be used only in RECV buffers\n");

    /* Before signaling the buffer is ready to receive data, destroy the data
     * header */
    packbuf_header_destroy(pb);

    switch (pb->transport) {
        case PB_GASPI: packbuf_gaspi_signal(pb); break;
        default: break;
    }
}
