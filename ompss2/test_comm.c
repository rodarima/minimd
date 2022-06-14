#include "test.h"

#define ENABLE_DEBUG 1
#include "log.h"
#include "types.h"

#include <string.h>

static int fake_natoms = 1;
static Vec fake_r = { 0.111111, 0.222222, 0.333333 };
static Vec fake_v = { 1.111111, 1.222222, 1.333333 };
static int fake_type = 123;

static void
check_header(PackBuf *pb, int magic, int iter)
{
    PackBufHeader *h = &pb->data->header;

    if (h->magic != magic)
        die("%s wrong magic %d\n", pb->name, h->magic);

    if (h->iter != iter)
        die("%s iter mismatch: recv %d, expected %d\n",
                pb->name, h->iter, iter);
}

static void
do_send(PackBuf *pb, enum pb_type type, enum pb_req req, int iter)
{
    if (type == PB_R) {
        packbuf_add(pb, &fake_r, NULL, NULL);
    } else if (type == PB_RT) {
        packbuf_add(pb, &fake_r, NULL, &fake_type);
    } else if (type == PB_RVT) {
        packbuf_add(pb, &fake_r, &fake_v, &fake_type);
    } else {
        die("unknown type\n");
    }

    pb->data->header.magic = PB_MAGIC_OK;
    pb->data->header.iter = iter;

	dbg("do_send: setting good magic in pb->data=%p\n", pb->data);

    check_header(pb, PB_MAGIC_OK, iter);

    dbg("sending %s iter=%d req=%s\n", pb->name, iter, PB_REQNAME(req));
    #pragma oss task
    packbuf_send(pb, req);
    #pragma oss taskwait

    check_header(pb, PB_MAGIC_OK, iter);
}

static void
do_recv(PackBuf *pb, enum pb_type type, enum pb_req req)
{
	PackBufHeader *h = &pb->data->header;

    packbuf_clear(pb);
    packbuf_grow(pb, 1);

    pb->natoms = fake_natoms;

	dbg("do_recv: setting bad magic in pb->data=%p\n", pb->data);

    //h->magic = PB_MAGIC_CLEAN;
    //h->iter = -1;
    //h->i = -1;

    //check_header(pb, PB_MAGIC_CLEAN, -1);

    dbg("receiving %s req=%s\n", pb->name, PB_REQNAME(req));

	dbg("do_recv: header BEFORE recv: magic=%d iter=%d icomm=%d\n",
		h->magic, h->iter, h->icomm);

    #pragma oss task
    packbuf_recv(pb, req);
    #pragma oss taskwait

    dbg("receiving %s req=%s has RELEASED dependencies\n",
		pb->name, PB_REQNAME(req));

	dbg("do_recv: header AFTER recv: magic=%d iter=%d icomm=%d\n",
		h->magic, h->iter, h->icomm);
}

static void
do_wait(PackBuf *pb, enum pb_type type, enum pb_req req)
{
    dbg("waiting %s req=%s\n", pb->name, PB_REQNAME(req));

	dbg("do_wait: testing magic for pb->data=%p dir=%s\n",
		pb->data, PB_DIRNAME(pb->dir));
    
    #pragma oss task
    packbuf_waitn(&pb, 1, req);

    #pragma oss taskwait
}

static void
do_check(PackBuf *pb, enum pb_type type, enum pb_req req, int iter)
{
    check_header(pb, PB_MAGIC_OK, iter);

    dbg("header seems OK: magic=%d iter=%d %s\n",
            pb->data->header.magic,
            pb->data->header.iter,
            pb->name);

    pb->data->header.magic = PB_MAGIC_CLEAN;
    pb->data->header.iter = -1;

    packbuf_clear(pb);
}

static void
do_barrier(int rank)
{
    fflush(stderr);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
        dbg("------------------------------------------------------------------------\n");

    MPI_Barrier(MPI_COMM_WORLD);
}

static void
clean_all_magic(Sim *sim, int iter)
{
    for (int ibox = 0; ibox < sim->nboxes; ibox++) {
        Box *box = &sim->box[ibox];
        for (int ineigh = 0; ineigh < NNEIGH; ineigh++) {
            Neigh *neigh = &box->neigh[ineigh];
            for (int type = 0; type < PB_NTYPES; type++) {
                for (int senddir = 0; senddir < PB_NDIR; senddir++) {
                    PackBuf *pb = &neigh->pb[type][senddir];
                    pb->data->header.magic = PB_MAGIC_CLEAN;
                    //pb->data->header.iter = -iter;
                    //pb->data->header.icomm = 1000 + sim->rank;
                }
            }
        }
    }
}

static void
clean_all_segments(Sim *sim)
{
    for (int ineigh = 0; ineigh < NNEIGH; ineigh++) {
        for (int idir = 0; idir < PB_NDIR; idir++) {
            void *ptr = sim->gaspi.buf[ineigh][idir];
            size_t len = sim->gaspi.segsize;

            memset(ptr, 0, len);
        }
    }
}

static void
find_nonzero(Sim *sim)
{
    for (int ineigh = 0; ineigh < NNEIGH; ineigh++) {
        for (int idir = 0; idir < PB_NDIR; idir++) {
            uint32_t *ptr = sim->gaspi.buf[ineigh][idir];
            size_t len = sim->gaspi.segsize;

            for (size_t i = 0; i < len; i++) {
                if (ptr[i] != 0)
                    dbg("found something at %p: %d\n", &ptr[i], ptr[i]);
            }

            memset(ptr, 0, len);
        }
    }
}

static void
search_magic(Sim *sim, PackBuf *sendpb, int iter)
{
    int ngood = 0;
    for (int ibox = 0; ibox < sim->nboxes; ibox++) {
        Box *box = &sim->box[ibox];
        for (int ineigh = 0; ineigh < NNEIGH; ineigh++) {
            Neigh *neigh = &box->neigh[ineigh];
            for (int type = 0; type < PB_NTYPES; type++) {
                for (int senddir = 0; senddir < PB_NDIR; senddir++) {
                    PackBuf *pb = &neigh->pb[type][senddir];
                    if (pb->data->header.magic == PB_MAGIC_OK) {
                        dbg("iter=%d found good magic in box=%d neigh=%d "
                                "type=%s dir=%s founditer=%d pb->data=%p\n",
                                iter, ibox, ineigh, PB_TYPENAME(type), PB_DIRNAME(senddir),
								pb->data->header.iter, pb->data);
                        ngood++;
                    }
                }
            }
        }
    }

    if (ngood == 0)
		dbg("ERR iter=%d rank=%d ngood = 0\n", iter, sim->rank);

    if (sendpb->transport != PB_SHM && ngood != 1)
        die("ngood=%d != 1 for internode\n", ngood);

    if (sendpb->transport == PB_SHM && ngood != 2)
        die("ngood=%d != 1 for intranode\n", ngood);
}

int
find_recv_endpoint(Sim *sim, Endpoint *src, Endpoint *dst)
{
    for (int ibox = 0; ibox < sim->nboxes; ibox++) {
        Box *box = &sim->box[ibox];
        for (int isenddir = 0; isenddir < NNEIGH; isenddir++) {
            Neigh *neigh = &box->neigh[isenddir];

            if (neigh->rank != src->rank)
                continue;

            if (neigh->boxid != src->ibox)
                continue;

            if (neigh->opposite->i != src->isenddir)
                continue;

            dst->rank = sim->rank;
            dst->ibox = ibox;
            dst->isenddir = isenddir;

            return 0;
        }
    }

    return -1;
}

/** Test one transaction at a time, iterating through each PackBuf. All
 * ranks wait until they have work to do. There is only one transaction
 * at any given moment. */
void
test_comm(Sim *sim)
{
    int iseq = 10000;

    for (int sendrank = 0; sendrank < sim->nranks; sendrank++) {
        int sending = 0;

        if (sim->rank == sendrank)
            sending = 1;

        for (int i = 0; i < sim->nboxes; i++) {
            Box *sendbox = &sim->box[i];

            for (int isendneigh = 0; isendneigh < NNEIGH; isendneigh++) {
                Neigh *sendneigh = &sendbox->neigh[isendneigh];

                for (int type = 0; type < PB_NTYPES; type++) {

                    PackBuf *sendpb = &sendneigh->pb[type][PB_SEND];

                    /* Ensure the PackBuf is the same in the box */
                    if (sendpb != sendbox->pb[type][PB_SEND][isendneigh])
                        die("box->pb=%p and neigh->pb=%p mismatch\n",
                                sendpb, sendbox->pb[type][PB_SEND][isendneigh]);

                    int receiving = 0;
                    Endpoint dst = { 0 };
                    Endpoint src = {
                        .rank = sendrank,
                        .ibox = i,
                        .isenddir = isendneigh
                    };

                    /* Compute the receiving endpoint */
                    if (find_recv_endpoint(sim, &src, &dst) == 0)
                        receiving = 1;

                    int recvrank = dst.rank;
                    Box *recvbox = &sim->box[dst.ibox];
                    Neigh *recvneigh = &recvbox->neigh[dst.isenddir];
                    PackBuf *recvpb = &recvneigh->pb[type][PB_RECV];

                    /* Ensure there is one and only one rank that
                     * performs the send and recv */

                    int nsenders = 0, nreceivers = 0;
                    MPI_Reduce(&sending, &nsenders, 1, MPI_INT,
                            MPI_SUM, 0, MPI_COMM_WORLD);
                    MPI_Reduce(&receiving, &nreceivers, 1, MPI_INT,
                            MPI_SUM, 0, MPI_COMM_WORLD);

                    if (sim->rank == 0) {
                        if (nsenders != 1)
                            die("wrong number of senders %d\n", nsenders);
                        if (nreceivers != 1)
                            die("wrong number of receivers %d\n", nreceivers);
                    }

                    if (sending && receiving && sendpb->transport != PB_SHM)
                        die("send transport should be SHM\n");

                    if (sending && receiving && recvpb->transport != PB_SHM)
                        die("recv transport should be SHM\n");

                    for (int req = 0; req < PB_NREQS; req++) {

                        if (req == PB_NATOMS
                                && sendpb->transport == PB_GASPI)
                            continue;

                        /* FIXME */
                        //if (sendpb->transport != PB_GASPI)
                        //    continue;

                        dbg("sendrank=%d sendbox=%d sendneigh=%d type=%s req=%s\n"
                                "  rank=%d req=%s iter=%d\n"
                                "  sending=%d (sendrank=%d sendbox=%d sendneigh=%2d curcoord=(%2d %2d %2d) delta=(%2d %2d %2d) dst coordw=(%2d %2d %2d))\n"
                                "  recving=%d (recvrank=%d recvbox=%d recvneigh=%2d curcoord=(%2d %2d %2d) delta=(%2d %2d %2d) dst coordw=(%2d %2d %2d))\n",
                                sendrank, sendbox->i, sendneigh->i, PB_TYPENAME(type), PB_REQNAME(req),
                                sim->rank, PB_REQNAME(req), iseq,
                                sending, sendrank, sendbox->i, sendneigh->i,
                                sendbox->idim[X],
                                sendbox->idim[Y],
                                sendbox->idim[Z],
                                sendneigh->delta[X],
                                sendneigh->delta[Y],
                                sendneigh->delta[Z],
                                sendneigh->boxcoordw[X],
                                sendneigh->boxcoordw[Y],
                                sendneigh->boxcoordw[Z],
                                receiving, recvrank, recvbox->i, recvneigh->i,
                                recvbox->idim[X],
                                recvbox->idim[Y],
                                recvbox->idim[Z],
                                recvneigh->delta[X],
                                recvneigh->delta[Y],
                                recvneigh->delta[Z],
                                recvneigh->boxcoordw[X],
                                recvneigh->boxcoordw[Y],
                                recvneigh->boxcoordw[Z]);

                        clean_all_magic(sim, iseq);

                        do_barrier(sim->rank);

                        if (sending)
                            do_send(sendpb, type, req, iseq);

                        if (receiving)
                            do_recv(recvpb, type, req);

                        if (sending)
                            do_wait(sendpb, type, req);

                        if (receiving)
                            do_wait(recvpb, type, req);

                        /* Check the magic of all PB, to see if there is
                         * another one which has been written */
                        if (receiving) {
                            search_magic(sim, sendpb, iseq);
                            //if (recvpb->transport == PB_GASPI)
                            //    find_nonzero(sim);
                        }

                        do_barrier(sim->rank);

                        if (receiving)
                            do_check(recvpb, type, req, iseq);

                        do_barrier(sim->rank);

                        if (sending)
                            do_check(sendpb, type, req, iseq);

                        do_barrier(sim->rank);

                        iseq++;
                    }
                }
            }
        }
    }
}
