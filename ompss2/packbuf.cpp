#include "types.h"

#include <stdlib.h>
#include <mpi.h>

#define PACKBUF_INCR 2000

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

//void
//packbuf_pack_sel_add(PackBuf *pb, Vec *src)
//{
//    for (int i = 0; i < pb->natoms; i++) {
//        for (int d = X; d <= Z; d++) {
//            pb->buf[i][d] = src[pb->sel[i]][d] + pb->add[d];
//        }
//    }
//}
//
//void
//packbuf_pack_r(PackBuf *pb, Vec *src)
//{
//    packbuf_grow(pb, 0);
//    packbuf_pack_sel_add(pb, src);
//}


void
packbuf_mpisend_buf(PackBuf *pb, int remoterank, int tag)
{
    fprintf(stderr, "packbuf sending %d atoms to rank %d\n",
            pb->natoms, remoterank);

    if (pb->natoms == 0)
        return;

    MPI_Isend((void *) pb->buf, pb->natoms * pb->atomsize, MPI_DOUBLE,
            remoterank, tag, pb->comm, &pb->req);
}

void
packbuf_mpisend(PackBuf *pb, int remoterank, int tag)
{
    fprintf(stderr, "packbuf sending %d atoms to rank %d\n",
            pb->natoms, remoterank);

    /* Send the number of atoms first */
    MPI_Send((void *) &pb->natoms, 1, MPI_INT,
            remoterank, tag, pb->comm);

    packbuf_mpisend_buf(pb, remoterank, tag);
}

void
packbuf_mpirecv_buf(PackBuf *pb, int remoterank, int tag, int natoms)
{
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
    if (src->natoms != 0) {
        packbuf_grow(dst, src->natoms);
        memcpy(dst->buf, src->buf, src->natoms * src->atomsize * sizeof(double));
    }
    dst->natoms = src->natoms;
}

/* FIXME: move to .h so the compiler can optimize the constant NULL
 * pointers */
void
packbuf_add(PackBuf *pb, Vec *r, Vec *v, int *type)
{
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

    if (j != pb->natoms * pb->atomsize) {
        fprintf(stderr, "packbuf_add atom size mismatch\n");
        abort();
    }
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
}

void
packbuf_unpack_sel(PackBuf *pb, Vec *r, Vec *v, int *types, int *sel)
{
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
}

void
packbuf_add_rt(PackBuf *pb, Vec r, int type)
{
    /* Ensure we have room for another atom */
    packbuf_grow_extra(pb, 1);

    int j = pb->natoms * pb->atomsize;

//    fprintf(stderr, "packing r %e %e %e into %p\n",
//            r[X], r[Y], r[Z], &pb->buf[j]);

    for (int d = X; d <= Z; d++)
        pb->buf[j++] = r[d];

    /* FIXME: We are sending the type as a double */
    pb->buf[j++] = (double) type;

    pb->natoms++;
}

void
packbuf_unpack_rt(PackBuf *pb, Vec *r, int *types)
{
    for (int i = 0, j = 0; i < pb->natoms; i++) {

        if (j+4 > pb->natoms * pb->atomsize)
            abort();

        for (int d = X; d <= Z; d++)
            r[i][d] = pb->buf[j++];

        types[i] = (int) pb->buf[j++];
    }
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
}
