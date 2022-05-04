#define ENABLE_DEBUG 0
#include "types.h"
#include "log.h"
#include "safe.h"
#include "packbuf.h"

#include <string.h>
#include <mpi.h>
#include <TAMPI.h>

static int
check_tag(int tag)
{
    /* Ensure the tag is within the MPI standard limit */
    if (tag >= 32767) {
        die("MPI tag exceed limit: %d >= %d\n", tag, 32767);
    }

    return tag;
}

void
packbuf_mpi_init(PackBuf *pb, int tag[PB_NREQS],
        int icomm, MPI_Comm *comm)
{
    pb->mode = PB_MPI;
    pb->mpi.icomm = icomm;
    pb->mpi.comm = comm;

    /* Allocate data with empty buf */
    pb->data = safe_calloc(1, sizeof(PackBufData));
    pb->nalloc = 0;

    for (int i = 0; i < PB_NREQS; i++) {
        pb->mpi.tag[i] = check_tag(tag[i]);
    }
}

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
    PackBufMPI *pbm = &pb->mpi;
    packbuf_switch(pb, PB_READY, PB_SENDING);

    int tag = pbm->tag[PB_BUF];

    dbg("send_buf: natoms=%d remoterank=%d tag=%d icomm=%d name='%s'\n",
            pb->natoms, pb->remoterank, tag, pbm->icomm,
            pb->name);

    if (ENABLE_NONBLOCKING_MPI) {
        if (pb->waitreq[PB_BUF])
            die("packbuf_mpi_send_buf: buffer in use\n");

        if (pb->natoms != 0) {

            isend((void *) pb->data->buf, pb->natoms * pb->atomsize, MPI_DOUBLE,
                    pb->remoterank, tag, *pbm->comm, &pbm->req[PB_BUF]);

            if (NEED_EXPLICIT_WAIT)
                pb->waitreq[PB_BUF] = 1;
        }
    } else {
        if (pb->natoms != 0) {
            MPI_Send((void *) pb->data->buf, pb->natoms * pb->atomsize,
                    MPI_DOUBLE, pb->remoterank, tag, *pbm->comm);
        }
    }

    packbuf_switch(pb, PB_SENDING, PB_READY);
}

static void
send_natoms(PackBuf *pb)
{
    int tag = pb->mpi.tag[PB_NATOMS];

    dbg("send_natoms: xnatoms=%d remoterank=%d tag=%d icomm=%d name='%s'\n",
            pb->natoms, pb->remoterank, tag, pb->mpi.icomm,
            pb->name);

    packbuf_switch(pb, PB_READY, PB_SENDING);

    pb->data->xnatoms = pb->natoms;

    void *buf = (void *) &pb->data->xnatoms;

    if (ENABLE_NONBLOCKING_MPI) {
        isend(buf, 1, MPI_INT, pb->remoterank, tag, *pb->mpi.comm,
                &pb->mpi.req[PB_NATOMS]);

        if (NEED_EXPLICIT_WAIT)
            pb->waitreq[PB_NATOMS] = 1;
    } else {
        MPI_Send(buf, 1, MPI_INT, pb->remoterank, tag, *pb->mpi.comm);
    }

    packbuf_switch(pb, PB_SENDING, PB_READY);
}

void
packbuf_mpi_send(PackBuf *pb, enum pb_req reqtype)
{
    if (pb->mode != PB_MPI)
        die("packbuf_mpi_send: incorrect mode\n");

    if (ENABLE_NONBLOCKING_MPI && pb->waitreq[reqtype]) {
        die("packbuf_mpi_send: buffer %s in use\n",
                PB_REQNAME(reqtype));
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
    int tag = pb->mpi.tag[PB_BUF];
    packbuf_switch(pb, PB_READY, PB_RECVING);

    dbg("recv_buf: recvnatoms=%d remoterank=%d tag=%d icomm=%d name='%s'\n",
            pb->data->xnatoms, pb->remoterank, tag, pb->mpi.icomm,
            pb->name);

    if (pb->mode != PB_MPI)
        die("packbuf_mpi_recv_buf: incorrect mode\n");

    if (ENABLE_NONBLOCKING_MPI && pb->waitreq[PB_BUF])
        die("packbuf_mpi_recv_buf: buffer in use\n");

    /* NOTE: Ensure xnatoms is properly set */
    int recvnatoms = pb->data->xnatoms;

    if (recvnatoms > 0) {
        /* Grow the buffer if needed */
        packbuf_grow(pb, recvnatoms);

        /* And receive that many atoms */
        int size = recvnatoms * pb->atomsize;

        if (ENABLE_NONBLOCKING_MPI) {
            irecv((void *) pb->data->buf, size, MPI_DOUBLE,
                    pb->remoterank, tag, *pb->mpi.comm, &pb->mpi.req[PB_BUF]);
            if (NEED_EXPLICIT_WAIT)
                pb->waitreq[PB_BUF] = 1;
        } else {
            MPI_Recv((void *) pb->data->buf, size, MPI_DOUBLE,
                    pb->remoterank, tag, *pb->mpi.comm, MPI_STATUS_IGNORE);
        }
    }

    packbuf_switch(pb, PB_RECVING, PB_READY);
}

/** Receive natoms into data.xnatoms. */
static void
recv_natoms(PackBuf *pb)
{
    int tag = pb->mpi.tag[PB_NATOMS];
    packbuf_switch(pb, PB_READY, PB_RECVING);
    /* Find out how many atoms I need to make room for */
    dbg("recv_natoms: natoms=? remoterank=%d tag=%d icomm=%d name='%s'\n",
            pb->remoterank, tag, pb->mpi.icomm, pb->name);

    if (pb->mode != PB_MPI)
        die("packbuf_mpi_recv_buf: incorrect mode\n");

    if (ENABLE_NONBLOCKING_MPI) {

        if (pb->waitreq[PB_NATOMS])
            die("packbuf_mpi_recv_natoms: buffer in use\n");

        irecv((void *) &pb->data->xnatoms, 1, MPI_INT,
                pb->remoterank, tag, *pb->mpi.comm, &pb->mpi.req[PB_NATOMS]);

        if (NEED_EXPLICIT_WAIT) {
            pb->waitreq[PB_NATOMS] = 1;
        }

    } else {
        MPI_Recv((void *) &pb->data->xnatoms, 1, MPI_INT,
                pb->remoterank, tag, *pb->mpi.comm, MPI_STATUS_IGNORE);
    }

    packbuf_switch(pb, PB_RECVING, PB_READY);
}

void
packbuf_mpi_recv(PackBuf *pb, enum pb_req reqtype)
{
    if (reqtype == PB_NATOMS) {
        recv_natoms(pb);
    } else {
        recv_buf(pb);
    }
}


#define MAXREQ NNEIGH

/**
 * Waits for all requests operating on pb->data->xnatoms to finish and copies the
 * result into pb->natoms.
 *
 * @param pbs The PackBuf array of pointers to wait from
 * @param n   Number of PackBuf elements in the array
 */
void
packbuf_mpi_waitn(PackBuf **pbs, int n, enum pb_req req)
{
    if (n > MAXREQ)
        die("too many requests");

    MPI_Request mpireq[MAXREQ];
    int nreq = 0;

    for (int i = 0; i < n; i++)
        packbuf_switch(pbs[i], PB_READY, PB_WAITING);

    if (ENABLE_SEQUENTIAL_MPIWAIT) {
        for (int i = 0; i < n; i++) {
            if (pbs[i]->waitreq[req])
                MPI_Wait(&pbs[i]->mpi.req[req], MPI_STATUS_IGNORE);
        }

    } else {

        for (int i = 0; i < n; i++) {
            if (pbs[i]->waitreq[req])
                memcpy(&mpireq[nreq++], &pbs[i]->mpi.req[req], sizeof(MPI_Request));
        }

        MPI_Waitall(nreq, mpireq, MPI_STATUSES_IGNORE);
    }

    for (int i = 0; i < n; i++) {
        pbs[i]->natoms = pbs[i]->data->xnatoms;
        pbs[i]->waitreq[req] = 0;
        packbuf_switch(pbs[i], PB_WAITING, PB_READY);
    }
}
