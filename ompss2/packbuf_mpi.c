#define ENABLE_DEBUG 1
#include "types.h"
#include "log.h"

#include <mpi.h>
#include <TAMPI.h>

static int
isend(const void *buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm, MPI_Request *request)
{
    if (ENABLE_NONBLOCKING_TAMPI) {
        return TAMPI_Isend(buf, count, datatype, dest, tag, comm, request);
    }

    return MPI_Isend(buf, count, datatype, dest, tag, comm, request);
}

static int
irecv(void *buf, int count, MPI_Datatype datatype, int source,
        int tag, MPI_Comm comm, MPI_Request *request)
{
    if (ENABLE_NONBLOCKING_TAMPI) {
        return TAMPI_Irecv(buf, count, datatype, source, tag, comm, request,
                MPI_STATUS_IGNORE);
    }

    return MPI_Irecv(buf, count, datatype, source, tag, comm, request);
}

void
packbuf_mpi_send_buf(PackBuf *pb)
{
    packbuf_switch(pb, PB_READY, PB_SENDING);

    dbg("packbuf_mpi_send_buf: natoms=%d remoterank=%d tag=%d icomm=%d\n",
            pb->natoms, pb->remoterank, pb->tag, pb->icomm);

    if (pb->gaspi)
        die("cannot use GASPI buffer with MPI\n");

    if (ENABLE_NONBLOCKING_MPI) {
        if (pb->waitreq)
            die("packbuf_mpi_send_buf: buffer in use\n");

        if (pb->natoms != 0) {

            isend((void *) pb->buf, pb->natoms * pb->atomsize, MPI_DOUBLE,
                    pb->remoterank, pb->tag, *pb->comm, &pb->req);

            if (NEED_EXPLICIT_WAIT)
                pb->waitreq = 1;
        }
    } else {
        if (pb->natoms != 0) {
            MPI_Send((void *) pb->buf, pb->natoms * pb->atomsize,
                    MPI_DOUBLE, pb->remoterank, pb->tag, *pb->comm);
        }
    }

    packbuf_switch(pb, PB_SENDING, PB_READY);
}

void
packbuf_mpi_send(PackBuf *pb)
{
    dbg("packbuf_mpi_send: natoms=%d remoterank=%d tag=%d icomm=%d\n",
            pb->natoms, pb->remoterank, pb->tag, pb->icomm);

    if (pb->gaspi)
        die("cannot use GASPI buffer with MPI\n");

    if (ENABLE_NONBLOCKING_MPI && pb->waitreqn)
        die("packbuf_mpi_send: buffer in use\n");

    void *buf = (void *) &pb->natoms;

    if (ENABLE_NONBLOCKING_MPI) {
        isend(buf, 1, MPI_INT, pb->remoterank, pb->tag, *pb->comm,
                &pb->reqn);

        if (NEED_EXPLICIT_WAIT)
            pb->waitreqn = 1;
    } else {
        MPI_Send(buf, 1, MPI_INT, pb->remoterank, pb->tag, *pb->comm);
    }

    packbuf_mpi_send_buf(pb);
}

void
packbuf_mpi_recv_buf(PackBuf *pb, int natoms)
{
    packbuf_switch(pb, PB_READY, PB_RECVING);

    dbg("packbuf_mpi_recv_buf: natoms=%d remoterank=%d tag=%d icomm=%d\n",
            natoms, pb->remoterank, pb->tag, pb->icomm);

    if (pb->gaspi)
        die("cannot use GASPI buffer with MPI\n");

    if (ENABLE_NONBLOCKING_MPI && pb->waitreq)
        die("packbuf_mpi_recv_buf: buffer in use\n");

    if (natoms > 0) {
        /* Grow the buffer if needed */
        packbuf_grow(pb, natoms);

        /* And receive that many atoms */
        int size = natoms * pb->atomsize;

        if (ENABLE_NONBLOCKING_MPI) {
            irecv((void *) pb->buf, size, MPI_DOUBLE,
                    pb->remoterank, pb->tag, *pb->comm, &pb->req);
            if (NEED_EXPLICIT_WAIT)
                pb->waitreq = 1;
        } else {
            MPI_Recv((void *) pb->buf, size, MPI_DOUBLE,
                    pb->remoterank, pb->tag, *pb->comm, MPI_STATUS_IGNORE);
        }
    }

    /* FIXME: this is dangerous as we are writing the natoms in the buffer
    while they may be still being written by MPI_Irecv */
    pb->natoms = natoms;
    packbuf_switch(pb, PB_RECVING, PB_READY);
}

void
packbuf_mpi_recv_natoms(PackBuf *pb)
{
    packbuf_switch(pb, PB_READY, PB_RECVING);
    /* Find out how many atoms I need to make room for */
    dbg("packbuf_mpi_recv_natoms: natoms=? remoterank=%d tag=%d icomm=%d\n",
            pb->remoterank, pb->tag, pb->icomm);

    if (pb->gaspi)
        die("cannot use GASPI buffer with MPI\n");

    if (ENABLE_NONBLOCKING_MPI) {

        if (pb->waitreqn)
            die("packbuf_mpi_recv_natoms: buffer in use\n");

        irecv((void *) &pb->recvnatoms, 1, MPI_INT,
                pb->remoterank, pb->tag, *pb->comm, &pb->reqn);

        if (NEED_EXPLICIT_WAIT) {
            pb->waitreqn = 1;
        }

    } else {
        MPI_Recv((void *) &pb->recvnatoms, 1, MPI_INT,
                pb->remoterank, pb->tag, *pb->comm, MPI_STATUS_IGNORE);
    }
    packbuf_switch(pb, PB_RECVING, PB_READY);
}
