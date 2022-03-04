#include "types.h"

#include <stdlib.h>
#include <mpi.h>

#define PACKBUF_INCR 2000

static void
packbuf_grow(PackBuf *pb, int n)
{
    if (n < pb->natoms)
        n = pb->natoms;

    if (pb->nalloc < n) {
        //if (pb->nalloc + PACKBUF_INCR < n)
        //    n = pb->nalloc + PACKBUF_INCR;

        size_t newalloc = sizeof(double) * pb->atomsize * n;
        //fprintf(stderr, "packbuf realloc to %ld bytes\n", newalloc);
        if (newalloc == 0)
            abort();
        pb->buf = (double *) safe_realloc(pb->buf, newalloc);
        pb->nalloc = n;
    }
}

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
packbuf_mpisend(PackBuf *pb, int remoterank, int tag)
{
    fprintf(stderr, "packbuf sending %d atoms to rank %d\n",
            pb->natoms, remoterank);

    /* Send the number of atoms first */
    MPI_Send((void *) &pb->natoms, 1, MPI_INT,
            remoterank, tag, pb->comm);

    if (pb->natoms == 0)
        return;

    MPI_Isend((void *) pb->buf, pb->natoms * pb->atomsize, MPI_DOUBLE,
            remoterank, tag, pb->comm, &pb->req);
}

void
packbuf_mpirecv(PackBuf *pb, int remoterank, int tag)
{
    /* Find out how many atoms I need to make room for */
    int natoms;
    MPI_Recv((void *) &natoms, 1, MPI_INT,
            remoterank, tag, pb->comm, MPI_STATUS_IGNORE);

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
packbuf_shmcopy(PackBuf *src, PackBuf *dst)
{
    if (src->natoms != 0) {
        packbuf_grow(dst, src->natoms);
        memcpy(src->buf, dst->buf, src->natoms * src->atomsize * sizeof(double));
    }
    dst->natoms = src->natoms;
}

void
packbuf_add_rvt(PackBuf *pb, Vec r, Vec v, int type)
{
    /* Ensure we have room for another atom */
    packbuf_grow_extra(pb, 1);

    int j = pb->natoms * pb->atomsize;

    for (int d = X; d <= Z; d++)
        pb->buf[j++] = r[d];

    for (int d = X; d <= Z; d++)
        pb->buf[j++] = v[d];

    /* FIXME: We are sending the type as a double */
    pb->buf[j++] = (double) type;

    pb->natoms++;
}

void
packbuf_unpack_rvt(PackBuf *pb, Vec *r, Vec *v, int *types)
{
    for (int i = 0, j = 0; i < pb->natoms; i++) {
        for (int d = X; d <= Z; d++)
            r[i][d] = pb->buf[j++];

        for (int d = X; d <= Z; d++)
            v[i][d] = pb->buf[j++];

        types[i] = (int) pb->buf[j++];
    }
}

void
packbuf_add_rt(PackBuf *pb, Vec r, int type)
{
    /* Ensure we have room for another atom */
    packbuf_grow_extra(pb, 1);

    int j = pb->natoms * pb->atomsize;

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
packbuf_init(PackBuf *pb, int atomsize)
{
    memset(pb, 0, sizeof(*pb));

    pb->atomsize = atomsize;
    MPI_Comm_dup(MPI_COMM_WORLD, &pb->comm);
}
