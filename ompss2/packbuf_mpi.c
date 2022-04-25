#define ENABLE_DEBUG 0
#include "types.h"
#include "log.h"
#include "packbuf.h"

#include <string.h>
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

static void
send_buf(PackBuf *pb)
{
    packbuf_switch(pb, PB_READY, PB_SENDING);

    int tag = pb->tag[PB_BUF];

    dbg("send_buf: natoms=%d remoterank=%d tag=%d icomm=%d name='%s'\n",
            pb->natoms, pb->remoterank, tag, pb->icomm,
            pb->name);

    if (ENABLE_NONBLOCKING_MPI) {
        if (pb->waitreq[PB_BUF])
            die("packbuf_mpi_send_buf: buffer in use\n");

        if (pb->natoms != 0) {

            isend((void *) pb->buf, pb->natoms * pb->atomsize, MPI_DOUBLE,
                    pb->remoterank, tag, *pb->comm, &pb->req[PB_BUF]);

            if (NEED_EXPLICIT_WAIT)
                pb->waitreq[PB_BUF] = 1;
        }
    } else {
        if (pb->natoms != 0) {
            MPI_Send((void *) pb->buf, pb->natoms * pb->atomsize,
                    MPI_DOUBLE, pb->remoterank, tag, *pb->comm);
        }
    }

    packbuf_switch(pb, PB_SENDING, PB_READY);
}

static void
send_natoms(PackBuf *pb)
{
    int tag = pb->tag[PB_NATOMS];

    dbg("send_natoms: natoms=%d remoterank=%d tag=%d icomm=%d name='%s'\n",
            pb->natoms, pb->remoterank, tag, pb->icomm,
            pb->name);

    packbuf_switch(pb, PB_READY, PB_SENDING);

    void *buf = (void *) &pb->natoms;

    if (ENABLE_NONBLOCKING_MPI) {
        isend(buf, 1, MPI_INT, pb->remoterank, tag, *pb->comm,
                &pb->req[PB_NATOMS]);

        if (NEED_EXPLICIT_WAIT)
            pb->waitreq[PB_NATOMS] = 1;
    } else {
        MPI_Send(buf, 1, MPI_INT, pb->remoterank, tag, *pb->comm);
    }

    packbuf_switch(pb, PB_SENDING, PB_READY);
}

void
packbuf_mpi_send(PackBuf *pb, enum pb_reqtype reqtype)
{
    if (pb->gaspi)
        die("cannot use GASPI buffer with MPI\n");

    if (ENABLE_NONBLOCKING_MPI && pb->waitreq[reqtype]) {
        die("packbuf_mpi_send: buffer %s in use\n",
                PB_REQTYPENAME(reqtype));
    }

    if (reqtype == PB_NATOMS) {
        send_natoms(pb);
    } else {
        send_buf(pb);
    }
}

static void
recv_buf(PackBuf *pb)
{
    int tag = pb->tag[PB_BUF];
    packbuf_switch(pb, PB_READY, PB_RECVING);

    dbg("recv_buf: recvnatoms=%d remoterank=%d tag=%d icomm=%d name='%s'\n",
            pb->recvnatoms, pb->remoterank, tag, pb->icomm,
            pb->name);

    if (pb->gaspi)
        die("cannot use GASPI buffer with MPI\n");

    if (ENABLE_NONBLOCKING_MPI && pb->waitreq[PB_BUF])
        die("packbuf_mpi_recv_buf: buffer in use\n");

    if (pb->recvnatoms > 0) {
        /* Grow the buffer if needed */
        packbuf_grow(pb, pb->recvnatoms);

        /* And receive that many atoms */
        int size = pb->recvnatoms * pb->atomsize;

        if (ENABLE_NONBLOCKING_MPI) {
            irecv((void *) pb->buf, size, MPI_DOUBLE,
                    pb->remoterank, tag, *pb->comm, &pb->req[PB_BUF]);
            if (NEED_EXPLICIT_WAIT)
                pb->waitreq[PB_BUF] = 1;
        } else {
            MPI_Recv((void *) pb->buf, size, MPI_DOUBLE,
                    pb->remoterank, tag, *pb->comm, MPI_STATUS_IGNORE);
        }
    }

    /* FIXME: this is dangerous as we are writing the natoms in the buffer
    while they may be still being written by MPI_Irecv. We can move this to the
    wait operation if we have NEED_EXPLICIT_WAIT. */
    pb->natoms = pb->recvnatoms;
    packbuf_switch(pb, PB_RECVING, PB_READY);
}

static void
recv_natoms(PackBuf *pb)
{
    int tag = pb->tag[PB_NATOMS];
    packbuf_switch(pb, PB_READY, PB_RECVING);
    /* Find out how many atoms I need to make room for */
    dbg("recv_natoms: natoms=? remoterank=%d tag=%d icomm=%d name='%s'\n",
            pb->remoterank, tag, pb->icomm, pb->name);

    if (pb->gaspi)
        die("cannot use GASPI buffer with MPI\n");

    if (ENABLE_NONBLOCKING_MPI) {

        if (pb->waitreq[PB_NATOMS])
            die("packbuf_mpi_recv_natoms: buffer in use\n");

        irecv((void *) &pb->recvnatoms, 1, MPI_INT,
                pb->remoterank, tag, *pb->comm, &pb->req[PB_NATOMS]);

        if (NEED_EXPLICIT_WAIT) {
            pb->waitreq[PB_NATOMS] = 1;
        }

    } else {
        MPI_Recv((void *) &pb->recvnatoms, 1, MPI_INT,
                pb->remoterank, tag, *pb->comm, MPI_STATUS_IGNORE);
    }

    packbuf_switch(pb, PB_RECVING, PB_READY);
}

void
packbuf_mpi_recv(PackBuf *pb, enum pb_reqtype reqtype)
{
    if (reqtype == PB_NATOMS) {
        recv_natoms(pb);
    } else {
        recv_buf(pb);
    }
}


#define MAXREQ NNEIGH

/**
 * Waits for all requests operating on "natoms" to finish by using MPI_Waitall
 *
 * @param pbs The PackBuf array of pointers to wait from
 * @param n   Number of PackBuf elements in the array
 */
void
packbuf_mpi_waitn(PackBuf **pbs, int n, enum pb_reqtype reqtype)
{
    if (n > MAXREQ)
        die("too many requests");

    MPI_Request req[MAXREQ];
    int nreq = 0;

    for (int i = 0; i < n; i++)
        packbuf_switch(pbs[i], PB_READY, PB_WAITING);

    if (ENABLE_SEQUENTIAL_MPIWAIT) {
        for (int i = 0; i < n; i++) {
            if (pbs[i]->waitreq[reqtype])
                MPI_Wait(&pbs[i]->req[reqtype], MPI_STATUS_IGNORE);
        }

    } else {

        for (int i = 0; i < n; i++) {
            if (pbs[i]->waitreq[reqtype])
                memcpy(&req[nreq++], &pbs[i]->req[reqtype], sizeof(MPI_Request));
        }

        MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE);
    }

    for (int i = 0; i < n; i++) {
        pbs[i]->waitreq[reqtype] = 0;
        packbuf_switch(pbs[i], PB_WAITING, PB_READY);
    }
}
