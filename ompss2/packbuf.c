#define ENABLE_DEBUG 0
#include "types.h"
#include "log.h"
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

/* Grows the buffer so that the allocated capacity can hold at least n
 * atoms. */
void
packbuf_grow(PackBuf *pb, int n)
{
    if (pb->gaspi)
        die("packbuf_grow cannot operate on GASPI buffer\n");

    if (n < pb->natoms)
        n = pb->natoms;

    if (pb->nalloc < n) {
        //if (pb->nalloc + PACKBUF_INCR < n)
        //    n = pb->nalloc + PACKBUF_INCR;

        size_t newalloc = sizeof(double) * pb->atomsize * n;

        if (newalloc == 0)
            abort();

        pb->buf = (double *) safe_realloc(pb->buf, newalloc);

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

void
packbuf_shmcopy(PackBuf *src, PackBuf *dst, enum pb_reqtype reqtype)
{
    dbg("packbuf_shmcopy: natoms=%d reqtype=%d\n",
            src->natoms, reqtype);

    packbuf_switch(src, PB_READY, PB_COPYING);
    packbuf_switch(dst, PB_READY, PB_COPYING);

    if (reqtype == PB_NATOMS)
        die("non-sense\n");

    if (src->natoms != 0) {
        packbuf_grow(dst, src->natoms);
        memcpy(dst->buf, src->buf, src->natoms * src->atomsize * sizeof(double));
    }

    dst->natoms = src->natoms;

    packbuf_switch(dst, PB_COPYING, PB_READY);
    packbuf_switch(src, PB_COPYING, PB_READY);
}

/* FIXME: move to .h so the compiler can optimize the constant NULL
 * pointers */
void
packbuf_add(PackBuf *pb, Vec *r, Vec *v, int *type)
{
    packbuf_switch(pb, PB_READY, PB_ADDING);
    /* Ensure we have room for another atom */
    packbuf_grow_extra(pb, 1);

    if (pb->waitreq[PB_BUF] || pb->waitreq[PB_NATOMS])
        die("packbuf_add: buffer in use\n");

    int j = pb->natoms * pb->atomsize;

    if (r != NULL) {
        for (int d = X; d <= Z; d++)
            pb->buf[j++] = r[0][d];
    }

    if (v != NULL) {
        for (int d = X; d <= Z; d++)
            pb->buf[j++] = v[0][d];
    }

    if (type != NULL) {
        /* FIXME: We are sending the type as a double */
        pb->buf[j++] = (double) type[0];
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

    if (pb->waitreq[PB_BUF] || pb->waitreq[PB_NATOMS])
        die("packbuf_unpack: buffer in use\n");

    for (int i = 0, j = 0; i < pb->natoms; i++) {
        if (r != NULL) {
            for (int d = X; d <= Z; d++)
                r[i][d] = pb->buf[j++];
        }

        if (v != NULL) {
            for (int d = X; d <= Z; d++)
                v[i][d] = pb->buf[j++];
        }

        if (types != NULL) {
            types[i] = (int) pb->buf[j++];
        }
    }
    packbuf_switch(pb, PB_UNPACKING, PB_READY);
}

void
packbuf_unpack_sel(PackBuf *pb, Vec *r, Vec *v, int *types, int *sel)
{
    packbuf_switch(pb, PB_READY, PB_UNPACKING);

    if (pb->waitreq[PB_BUF] || pb->waitreq[PB_NATOMS])
        die("packbuf_unpack_sel: buffer in use\n");

    for (int i = 0, j = 0; i < pb->natoms; i++) {
        if (r != NULL) {
            for (int d = X; d <= Z; d++)
                r[sel[i]][d] = pb->buf[j++];
        }

        if (v != NULL) {
            for (int d = X; d <= Z; d++)
                v[sel[i]][d] = pb->buf[j++];
        }

        if (types != NULL) {
            types[sel[i]] = (int) pb->buf[j++];
        }
    }
    packbuf_switch(pb, PB_UNPACKING, PB_READY);
}

void
packbuf_clear(PackBuf *pb)
{
    packbuf_switch(pb, PB_READY, PB_CLEANING);

    if (pb->waitreq[PB_BUF] || pb->waitreq[PB_NATOMS])
        die("packbuf_unpack_sel: buffer in use\n");

    pb->natoms = 0;

    packbuf_switch(pb, PB_CLEANING, PB_READY);
}

void
packbuf_send(PackBuf *pb, enum pb_reqtype reqtype)
{
    if (ENABLE_GASPI && pb->gaspi) {
        packbuf_gaspi_send(pb, reqtype);
    } else {
        packbuf_mpi_send(pb, reqtype);
    }
}

void
packbuf_recv(PackBuf *pb, enum pb_reqtype reqtype)
{
    if (ENABLE_GASPI && pb->gaspi) {
        packbuf_gaspi_recv(pb, reqtype);
    } else {
        packbuf_mpi_recv(pb, reqtype);
    }
}

void
packbuf_init(PackBuf *pb, int enable_sel, int atomsize,
        int remoterank, int tag[PB_NREQTYPES], int icomm, MPI_Comm *comm)
{
    memset(pb, 0, sizeof(*pb));

    pb->atomsize = atomsize;
    pb->enable_sel = enable_sel;
    pb->icomm = icomm;
    pb->comm = comm;

    if (remoterank < 0)
        abort();

    pb->remoterank = remoterank;

    for (int i = 0; i < PB_NREQTYPES; i++)
        pb->tag[i] = tag[i];

    packbuf_switch(pb, PB_GARBAGE, PB_READY);
}

enum pb_type
packbuf_opposite_dir(enum pb_type type)
{
    switch (type) {
        case PB_SEND_R: return PB_RECV_R;
        case PB_SEND_RT: return PB_RECV_RT;
        case PB_SEND_RVT: return PB_RECV_RVT;
        case PB_RECV_R: return PB_SEND_R;
        case PB_RECV_RT: return PB_SEND_RT;
        case PB_RECV_RVT: return PB_SEND_RVT;
        default: die("unknown pb_type\n");
    }

    /* Not reached */
    return 0;
}
