#define ENABLE_DEBUG 0
#include "types.h"
#include "log.h"
#include "safe.h"
#include "packbuf.h"

#include <string.h>
#include <stdlib.h>

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
packbuf_send(PackBuf *pb, enum pb_req reqtype)
{
    if (pb->in_transfer[reqtype])
        die("packbuf_send: already transferring data %s\n", pb->name);

    switch (pb->transport) {
        case PB_GASPI: packbuf_gaspi_send(pb, reqtype); break;
        case PB_MPI: packbuf_mpi_send(pb, reqtype); break;
        case PB_SHM: packbuf_shm_send(pb, reqtype); break;
        default: die("packbuf_send: bad transport\n");
    }
}

static void
destroy_header(PackBuf *pb)
{
    PackBufHeader *h = &pb->data->header;

    if (pb->in_transfer[PB_BUF] || pb->in_transfer[PB_NATOMS])
        die("cannot destroy header: in transfer is set %s", pb->name);

    if (h->magic != PB_MAGIC_OK) {
        die("cannot destroy header: magic is already clean %s\n",
                pb->name);
    }

    memset(h, 0xff, sizeof(*h));

    h->magic = PB_MAGIC_CLEAN;

    if (pb->nalloc >= 1 && pb->atomsize >= 1) {
        pb->data->buf[0] = 0.1234567;
    }
}

void
packbuf_recv(PackBuf *pb, enum pb_req reqtype)
{
    if (pb->in_transfer[reqtype])
        die("packbuf_recv: already transferring data %s\n", pb->name);

    /* Before receiving data, destroy the data header */
    //destroy_header(pb);

    switch (pb->transport) {
        case PB_GASPI: packbuf_gaspi_recv(pb, reqtype); break;
        case PB_MPI: packbuf_mpi_recv(pb, reqtype); break;
        case PB_SHM: packbuf_shm_recv(pb, reqtype); break;
        default: die("packbuf_recv: bad transport\n");
    }
}

void
packbuf_init(PackBuf *pb, enum pb_dir dir, int ineigh, int enable_sel, int atomsize, int remoterank)
{
    memset(pb, 0, sizeof(*pb));

    pb->atomsize = atomsize;
    pb->enable_sel = enable_sel;
    pb->transport = PB_BAD;
    pb->dir = dir;
    pb->ineigh = ineigh;

    if (remoterank < 0)
        die("packbuf_init: negative remote rank %d\n", remoterank);

    pb->remoterank = remoterank;
    pb->data = NULL;

    packbuf_switch(pb, PB_GARBAGE, PB_READY);
}

/** Ensure the header matches with the expected values */
void
packbuf_check_header(PackBuf *pb, int iter)
{
    PackBufHeader *h = &pb->data->header;

    if (h->magic != PB_MAGIC_OK)
        die("%s wrong magic %d\n", pb->name, h->magic);

    if (iter != -666 && h->iter != iter)
        die("%s iter mismatch: recv %d, expected %d\n",
                pb->name, h->iter, iter);

    if (h->dstbox != pb->box)
        die("%s box mismatch: recv %d, expected %d\n",
                pb->name, h->dstbox, pb->box);

    if (h->senddir != pb->senddir)
        die("%s senddir mismatch: recv %d, expected %d\n",
                pb->name, h->senddir, pb->senddir);

    if (pb->transport == PB_MPI) {
        if (h->icomm != pb->mpi.icomm)
            die("%s icomm mismatch: recv %d, expected %d\n",
                    pb->name, h->icomm, pb->mpi.icomm);
    }

    dbg("%s header ok\n", pb->name);
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
