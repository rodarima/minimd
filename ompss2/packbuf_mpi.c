#define ENABLE_DEBUG 1
#include "types.h"
#include "log.h"
#include "safe.h"
#include "packbuf.h"
#include "trace.h"

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
    pb->transport = PB_MPI;
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
isend(PackBuf *pb, const void *buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm, MPI_Request *request)
{
    char label[1024];
    PackBufHeader *h = &pb->data->header;
    sprintf(label, "isend xnatoms=%d count=%d dest=%d tag=%d comm=%d (%s)\n"
        "header: magic=%d iter=%d srcbox=%d senddir=%d dstbox=%d icomm=%d",
        pb->data->xnatoms, count, dest, tag, comm, pb->name,
        h->magic, h->iter, h->srcbox, h->senddir, h->dstbox, h->icomm);
    trace_event(pb->box * NNEIGH + pb->ineigh, label, "#00ffff");

    if (ENABLE_NONBLOCKING_TAMPI) {
        return TAMPI_Isend(buf, count, datatype, dest, tag, comm, request);
    }

    return MPI_Isend(buf, count, datatype, dest, tag, comm, request);
}

static int
irecv(PackBuf *pb, void *buf, int count, MPI_Datatype datatype, int source,
        int tag, MPI_Comm comm, MPI_Request *request)
{
    char label[1024];
    PackBufHeader *h = &pb->data->header;
    sprintf(label, "irecv count=%d source=%d tag=%d comm=%d (%s)\n"
        "header: magic=%d iter=%d srcbox=%d senddir=%d dstbox=%d icomm=%d",
        count, source, tag, comm, pb->name,
        h->magic, h->iter, h->srcbox, h->senddir, h->dstbox, h->icomm);
    trace_event(pb->box * NNEIGH + pb->ineigh, label, "#00ff00");

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

    /* Skip the double buf[] at the end */
    int bytes = sizeof(*pb->data) +
        pb->natoms * pb->atomsize * sizeof(double);
    void *buf = pb->data;

    char label[1024];
    sprintf(label, "packbuf_mpi.c:send_buf natoms=%d %s", pb->natoms, pb->name);
    trace_event(pb->box * NNEIGH + pb->ineigh, label, "#0000ff");

    if (ENABLE_NONBLOCKING_MPI) {
        if (pb->mpi.waitreq[PB_BUF])
            die("packbuf_mpi_send_buf: buffer in use\n");

        if (pb->natoms != 0) {

            isend(pb, buf, bytes, MPI_BYTE,
                    pb->remoterank, tag, *pbm->comm, &pbm->req[PB_BUF]);

            if (NEED_EXPLICIT_WAIT)
                pb->mpi.waitreq[PB_BUF] = 1;
        }
    } else {
        if (pb->natoms != 0) {
            MPI_Send(buf, bytes, MPI_BYTE,
                    pb->remoterank, tag, *pbm->comm);
        }
    }

    if (pb->natoms != 0)
        pb->in_transfer[PB_BUF] = 1;

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

    /* Skip the double buf[] at the end */
    int bytes = sizeof(*pb->data);
    void *buf = pb->data;

    if (ENABLE_NONBLOCKING_MPI) {
        isend(pb, buf, bytes, MPI_BYTE,
                pb->remoterank, tag, *pb->mpi.comm,
                &pb->mpi.req[PB_NATOMS]);

        if (NEED_EXPLICIT_WAIT)
            pb->mpi.waitreq[PB_NATOMS] = 1;
    } else {
        MPI_Send(buf, bytes, MPI_BYTE,
                pb->remoterank, tag, *pb->mpi.comm);
    }

    pb->in_transfer[PB_NATOMS] = 1;

    packbuf_switch(pb, PB_SENDING, PB_READY);
}

void
packbuf_mpi_send(PackBuf *pb, enum pb_req reqtype)
{
    if (pb->transport != PB_MPI)
        die("packbuf_mpi_send: incorrect transport\n");

    if (pb->in_transfer[reqtype])
        die("packbuf_mpi_send: already transferring data %s\n", pb->name);

    if (ENABLE_NONBLOCKING_MPI && pb->mpi.waitreq[reqtype]) {
        die("packbuf_mpi_send: buffer %s in use\n",
                PB_REQNAME(reqtype));
    }

    pb->data->xnatoms = pb->natoms;

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

    if (pb->transport != PB_MPI)
        die("packbuf_mpi_recv_buf: incorrect transport\n");

    /* Use pb->natoms as the number of atoms to be received */
    int recvnatoms = pb->natoms;

    dbg("recv_buf: recvnatoms=%d remoterank=%d tag=%d icomm=%d name='%s'\n",
            pb->data->xnatoms, pb->remoterank, tag, pb->mpi.icomm,
            pb->name);

    if (ENABLE_NONBLOCKING_MPI && pb->mpi.waitreq[PB_BUF])
        die("packbuf_mpi_recv_buf: buffer in use\n");

    if (recvnatoms > 0) {
        /* Grow the buffer if needed */
        packbuf_grow(pb, recvnatoms);

        /* And receive the header with the atom data */
        int bytes = sizeof(*pb->data)
            + recvnatoms * pb->atomsize * sizeof(double);

        void *buf = pb->data;

        if (ENABLE_NONBLOCKING_MPI) {
            irecv(pb, buf, bytes, MPI_BYTE,
                    pb->remoterank, tag, *pb->mpi.comm,
                    &pb->mpi.req[PB_BUF]);

            if (NEED_EXPLICIT_WAIT)
                pb->mpi.waitreq[PB_BUF] = 1;
        } else {
            MPI_Recv(buf, bytes, MPI_BYTE,
                    pb->remoterank, tag, *pb->mpi.comm,
                    MPI_STATUS_IGNORE);
        }

        pb->in_transfer[PB_BUF] = 1;
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

    if (pb->transport != PB_MPI)
        die("packbuf_mpi_recv_buf: incorrect transport\n");

    /* And receive the header with the atom data */
    int bytes = sizeof(*pb->data);
    void *buf = pb->data;

    if (ENABLE_NONBLOCKING_MPI) {

        if (pb->mpi.waitreq[PB_NATOMS])
            die("packbuf_mpi_recv_natoms: buffer in use\n");

        irecv(pb, buf, bytes, MPI_BYTE,
                pb->remoterank, tag, *pb->mpi.comm,
                &pb->mpi.req[PB_NATOMS]);

        if (NEED_EXPLICIT_WAIT)
            pb->mpi.waitreq[PB_NATOMS] = 1;

    } else {
        MPI_Recv(buf, bytes, MPI_BYTE,
                pb->remoterank, tag, *pb->mpi.comm,
                MPI_STATUS_IGNORE);

        dbg("recv_natoms: MPI_Recv natoms=%d %s\n",
                pb->data->xnatoms, pb->name);
    }

    pb->in_transfer[PB_NATOMS] = 1;

    packbuf_switch(pb, PB_RECVING, PB_READY);
}

void
packbuf_mpi_recv(PackBuf *pb, enum pb_req reqtype)
{
    if (pb->transport != PB_MPI)
        die("packbuf_mpi_recv: incorrect transport %s\n", pb->name);

    if (pb->in_transfer[reqtype])
        die("packbuf_mpi_recv: already transferring data %s\n", pb->name);

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
    MPI_Request mpireq[MAXREQ];
    int nreq = 0;

    for (int i = 0; i < n; i++) {
        PackBuf *pb = pbs[i];

        /* Ignore PackBufs which are not in transfer (natoms=0) */
        if (pb->transport != PB_MPI || !pb->in_transfer[req])
            continue;

        packbuf_switch(pb, PB_READY, PB_WAITING);

        /* No need to do any wait */
        if (!NEED_EXPLICIT_WAIT)
            continue;

        if (pb->in_transfer[req] && !pb->mpi.waitreq[req])
            die("packbuf_mpi_waitn: waitreq is not set %s\n", pb->name);

        if (nreq >= MAXREQ)
            die("packbuf_mpi_waitn: too many requests");

        if (ENABLE_SEQUENTIAL_MPIWAIT)
            MPI_Wait(&pb->mpi.req[req], MPI_STATUS_IGNORE);
        else
            memcpy(&mpireq[nreq++], &pb->mpi.req[req], sizeof(MPI_Request));
    }

    if (NEED_EXPLICIT_WAIT && !ENABLE_SEQUENTIAL_MPIWAIT && nreq != 0)
        MPI_Waitall(nreq, mpireq, MPI_STATUSES_IGNORE);

    for (int i = 0; i < n; i++) {
        PackBuf *pb = pbs[i];

        if (pb->transport != PB_MPI || !pb->in_transfer[req])
            continue;

        if (pb->dir == PB_RECV) {
            packbuf_check_header(pb, -666);
            pb->natoms = pb->data->xnatoms;
        }

        /* Clear all in_transfer and waitreq flags */
        pb->mpi.waitreq[req] = 0;
        pb->in_transfer[req] = 0;

        packbuf_switch(pb, PB_WAITING, PB_READY);
    }
}
