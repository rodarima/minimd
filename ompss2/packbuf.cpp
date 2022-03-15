#include "types.h"

#include <stdlib.h>
#include <mpi.h>

#define PACKBUF_INCR 2000

static void
packbuf_switch(PackBuf *pb, enum packbuf_state prev, enum packbuf_state next)
{
    if (pb->state != prev) {
        fprintf(stderr, "error: packbuf in state %d, expected %d\n",
                pb->state, prev);
        abort();
    }

    pb->state = next;
}

/* Grows the buffer so that the allocated capacity can hold at least n
 * atoms. */
static void
packbuf_grow(PackBuf *pb, int n)
{
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
packbuf_mpisend_buf(PackBuf *pb, int remoterank, int tag)
{
    packbuf_switch(pb, PB_READY, PB_SENDING);

    fprintf(stderr, "packbuf sending %d atoms to rank %d\n",
            pb->natoms, remoterank);

    if (remoterank < 0)
        abort();

    if (pb->natoms == 0)
        return;

    MPI_Isend((void *) pb->buf, pb->natoms * pb->atomsize, MPI_DOUBLE,
            remoterank, tag, pb->comm, &pb->req);

    /* FIXME: Wait for the communication to finish */

    packbuf_switch(pb, PB_SENDING, PB_READY);
}

void
packbuf_mpisend(PackBuf *pb, int remoterank, int tag)
{
    fprintf(stderr, "packbuf sending %d atoms to rank %d\n",
            pb->natoms, remoterank);

    if (remoterank < 0)
        abort();

    /* Send the number of atoms first */
    MPI_Send((void *) &pb->natoms, 1, MPI_INT,
            remoterank, tag, pb->comm);

    packbuf_mpisend_buf(pb, remoterank, tag);
}

void
packbuf_mpirecv_buf(PackBuf *pb, int remoterank, int tag, int natoms)
{
    packbuf_switch(pb, PB_READY, PB_RECVING);
    fprintf(stderr, "packbuf receiving %d atoms from rank %d\n",
            natoms, remoterank);

    if (natoms > 0) {
        /* Grow the buffer if needed */
        packbuf_grow(pb, natoms);

        /* And received that many atoms */
        int size = natoms * pb->atomsize;

        //pb->ready = 0;
        MPI_Irecv((void *) pb->buf, size, MPI_DOUBLE,
                remoterank, tag, pb->comm, &pb->req);
    }

    pb->natoms = natoms;
    packbuf_switch(pb, PB_RECVING, PB_READY);
}

void
packbuf_mpirecv(PackBuf *pb, int remoterank, int tag)
{
    /* Find out how many atoms I need to make room for */
    int natoms;
    MPI_Recv((void *) &natoms, 1, MPI_INT,
            remoterank, tag, pb->comm, MPI_STATUS_IGNORE);

    packbuf_mpirecv_buf(pb, remoterank, tag, natoms);
}

void
packbuf_shmcopy(PackBuf *src, PackBuf *dst)
{
    packbuf_switch(src, PB_READY, PB_COPYING);
    packbuf_switch(dst, PB_READY, PB_COPYING);

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

//    fprintf(stderr, "packbuf %p now has %d atoms (alloc %d)\n",
//            pb, pb->natoms, pb->nalloc);

    if (j != pb->natoms * pb->atomsize) {
        fprintf(stderr, "packbuf_add atom size mismatch\n");
        abort();
    }

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
    pb->natoms = 0;
}

void
packbuf_init(PackBuf *pb, int enable_sel, int atomsize)
{
    memset(pb, 0, sizeof(*pb));

    pb->atomsize = atomsize;
    pb->enable_sel = enable_sel;
    MPI_Comm_dup(MPI_COMM_WORLD, &pb->comm);
    packbuf_switch(pb, PB_GARBAGE, PB_READY);
}
