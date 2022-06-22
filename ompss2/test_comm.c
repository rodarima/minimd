#include "test.h"

#define ENABLE_DEBUG 0
#include "log.h"
#include "types.h"

#include <string.h>

static int fake_natoms = 1;
static Vec fake_r = { 0.111111, 0.222222, 0.333333 };
static Vec fake_v = { 1.111111, 1.222222, 1.333333 };
static int fake_type = 123;

static void
check_header(PackBuf *pb, int magic)
{
    PackBufHeader *h = &pb->data->header;

    if (h->magic != magic)
        die("Wrong magic %d: %s\n", h->magic, pb->name);
}

static void
do_send(PackBuf *pb, enum pb_req req)
{
    if (pb->type == PB_R) {
        packbuf_add(pb, &fake_r, NULL, NULL);
    } else if (pb->type == PB_RT) {
        packbuf_add(pb, &fake_r, NULL, &fake_type);
    } else if (pb->type == PB_RVT) {
        packbuf_add(pb, &fake_r, &fake_v, &fake_type);
    } else {
        die("unknown type\n");
    }

    dbg("sending %s iseq=%d req=%s\n", pb->name,
            pb->data->header.iseq, PB_REQNAME(req));

    #pragma oss task
    packbuf_send(pb, req);
    #pragma oss taskwait

    check_header(pb, PB_MAGIC_OK);
}

static void
do_recv(PackBuf *pb, enum pb_req req)
{
    PackBufHeader *h = &pb->data->header;

    packbuf_clear(pb);
    packbuf_grow(pb, 1);

    pb->natoms = fake_natoms;

    dbg("receiving %s req=%s\n", pb->name, PB_REQNAME(req));

    dbg("do_recv: header BEFORE recv: magic=%d iseq=%d\n",
        h->magic, h->iseq);

    #pragma oss task
    packbuf_recv(pb, req);
    #pragma oss taskwait

    dbg("receiving %s req=%s has RELEASED dependencies\n",
        pb->name, PB_REQNAME(req));

    dbg("do_recv: header AFTER recv: magic=%d iseq=%d\n",
        h->magic, h->iseq);

    /* The data may still not be available yet, until wait finishes */
}

static void
do_wait(PackBuf *pb, enum pb_req req)
{
    dbg("waiting %s req=%s\n", pb->name, PB_REQNAME(req));

    dbg("do_wait: testing magic for pb->data=%p dir=%s\n",
        pb->data, PB_DIRNAME(pb->dir));
    
    #pragma oss task
    packbuf_waitn(&pb, 1, req);

    #pragma oss taskwait
}

static void
do_check(PackBuf *pb, enum pb_req req)
{
    check_header(pb, PB_MAGIC_OK);

    dbg("header seems OK: magic=%d iseq=%d %s\n",
            pb->data->header.magic,
            pb->data->header.iseq,
            pb->name);

    packbuf_header_destroy(pb);

    dbg("ensure the header got destroyed\n");
    check_header(pb, PB_MAGIC_DESTROYED);
    dbg("header seems destroyed\n");

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
clean_all_magic(Sim *sim, int magic)
{
    for (int ibox = 0; ibox < sim->nboxes; ibox++) {
        Box *box = &sim->box[ibox];
        for (int ineigh = 0; ineigh < NNEIGH; ineigh++) {
            Neigh *neigh = &box->neigh[ineigh];
            for (int type = 0; type < PB_NTYPES; type++) {
                for (int senddir = 0; senddir < PB_NDIR; senddir++) {
                    PackBuf *pb = &neigh->pb[type][senddir];
                    pb->data->header.magic = magic;
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
search_magic(Sim *sim, int use_shm)
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
                        dbg("found good magic in box=%d neigh=%d "
                                "type=%s dir=%s iseq=%d pb->data=%p\n",
                                ibox, ineigh,
                                PB_TYPENAME(type), PB_DIRNAME(senddir),
                                pb->data->header.iseq, pb->data);
                        ngood++;
                    }
                }
            }
        }
    }

    if (ngood == 0)
        dbg("ERR rank=%d ngood = 0\n", sim->rank);

    int nexpected = use_shm ? 2 : 1;
    if (ngood != nexpected) {
        die("ngood=%d != %d with use_shm=%d\n",
                ngood, nexpected, use_shm);
    }
}

int
find_recv_endpoint(Sim *sim, Endpoint *src, Endpoint *dst)
{
    for (int ibox = 0; ibox < sim->nboxes; ibox++) {
        Box *box = &sim->box[ibox];

        Neigh *sendneigh = box->neigh[src->index].opposite;

        if (sendneigh->rank != src->rank)
            continue;

        if (sendneigh->boxid != src->box)
            continue;

        dst->rank = sim->rank;
        dst->box = ibox;
        dst->index = src->index;
        dst->type = src->type;
        dst->dir = PB_RECV;
        endpoint_set_name(dst);

        dbg("found matching endpoint %s -> %s\n", src->name, dst->name);

        return 0;
    }

    return -1;
}

/** Ensure there is one and only one rank that performs the send and recv */
static void
check_uniqueness(Sim *sim, int sending, int receiving)
{
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
}

static void
do_transactions(Sim *sim, PackBuf *srcpb, PackBuf *dstpb)
{
    int sending = (srcpb != NULL);
    int receiving = (dstpb != NULL);

    for (int req = 0; req < PB_NREQS; req++) {

        /* Don't test GASPI with NATOMS, as it's not implemented
         * (no need) */
        if (sending && req == PB_NATOMS && srcpb->transport == PB_GASPI)
            continue;

        if (receiving && req == PB_NATOMS && dstpb->transport == PB_GASPI)
            continue;

        if (sending)
            check_header(srcpb, PB_MAGIC_DESTROYED);

        if (receiving)
            check_header(dstpb, PB_MAGIC_DESTROYED);

        do_barrier(sim->rank);

        if (sending)
            do_send(srcpb, req);

        if (receiving)
            do_recv(dstpb, req);

        if (sending)
            do_wait(srcpb, req);

        if (receiving)
            do_wait(dstpb, req);

        int use_shm = (sending && receiving);

        /* Check the magic of all PB, to see if there is
         * another one which has been written */
        if (receiving)
            search_magic(sim, use_shm);

        do_barrier(sim->rank);

        if (receiving)
            do_check(dstpb, req);

        do_barrier(sim->rank);

        if (sending)
            do_check(srcpb, req);

        do_barrier(sim->rank);
    }
}

/** Test one transaction at a time, iterating through each PackBuf. All
 * ranks wait until they have work to do. There is only one transaction
 * at any given moment. */
static void
test_comm_run(Sim *sim, int run)
{
    /* First we build the source endpoint, and from that we compute the
     * destination endpoint. If the current rank matches any of those
     * endpoints then it perform the send, receive or both transactions
     * (shared memory) */

    Endpoint src = { 0 };

    src.dir = PB_SEND;

    for (src.rank = 0; src.rank < sim->nranks; src.rank++) {
        int sending = 0;

        if (src.rank == sim->rank)
            sending = 1;

        for (src.box = 0; src.box < sim->nboxes; src.box++) {
            Box *srcbox = &sim->box[src.box];
            err("test_comm_run run=%d src.rank=%d src.box=%d\n",
                    run, src.rank, src.box);

            for (src.index = 0; src.index < NNEIGH; src.index++) {
                Neigh *srcneigh = &srcbox->neigh[src.index];

                for (src.type = 0; src.type < PB_NTYPES; src.type++) {
                    PackBuf *srcpb = NULL;
                    PackBuf *dstpb = NULL;

                    endpoint_set_name(&src);

                    if (sending) {

                        srcpb = &srcneigh->pb[src.type][src.dir];

                        /* Ensure the PackBuf is the same in the box */
                        if (srcpb != srcbox->pb[src.type][src.dir][src.index]) {
                            die("box->pb=%p and neigh->pb=%p mismatch\n",
                                    srcpb, srcbox->pb[src.type][src.dir][src.index]);
                        }
                    }

                    int receiving = 0;
                    Endpoint dst = { 0 };

                    /* Compute the receiving endpoint */
                    if (find_recv_endpoint(sim, &src, &dst) == 0)
                        receiving = 1;

                    if (receiving) {
                        dstpb = endpoint_get_packbuf(sim, &dst);

                        if (dstpb == NULL)
                            die("bad dstpb\n");

                        /* Ensure the remote endpoint of the PackBuf
                         * matches with the source endpoint we are
                         * testing */

                        if (!endpoint_is_same(&src, &dstpb->remote)) {
                            die("Wrong PackBuf remote endpoint:\n"
                                    "  found:    %s\n"
                                    "  expected: %s\n",
                                    dstpb->remote.name, src.name);
                        }
                    }

                    check_uniqueness(sim, sending, receiving);

                    if (sending && receiving && srcpb->transport != PB_SHM)
                        die("send transport should be SHM\n");

                    if (sending && receiving && dstpb->transport != PB_SHM)
                        die("dst transport should be SHM\n");

                    do_transactions(sim, srcpb, dstpb);
                }
            }
        }
    }
}

void
test_comm(Sim *sim)
{
    /* Test thoroughly all communications several times */
    int nruns = 3;

    for (int run = 0; run < nruns; run++) {
        test_comm_run(sim, run);
    }
}
