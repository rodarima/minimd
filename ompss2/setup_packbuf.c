#define ENABLE_DEBUG 0
#include "log.h"
#include "types.h"
#include "gaspi_check.h"

#include <GASPI.h>
#include <TAGASPI.h>
#include <stdlib.h>

static int
get_isegment(int ineigh, int dir)
{
    return ineigh * PB_NDIR + dir;
}

static void
setup_gaspi_segment(Sim *sim, int ineigh, int dir, size_t nbytes)
{
    void *seg;

    if ((seg = malloc(nbytes)) == NULL)
        die("malloc of %zu bytes failed\n", nbytes);

    sim->gaspi.buf[ineigh][dir] = seg;

    int iseg = get_isegment(ineigh, dir);

    CHECK(gaspi_segment_use(iseg, seg, nbytes,
                GASPI_GROUP_ALL, GASPI_BLOCK, 0));

    CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
}

static void
setup_gaspi_segments(Sim *sim)
{
    /* Setup GASPI config */
    gaspi_config_t conf;
    CHECK(gaspi_config_get(&conf));
    conf.build_infrastructure = GASPI_TOPOLOGY_DYNAMIC;
    conf.queue_size_max = 4*1024;
    CHECK(gaspi_config_set(conf));

    CHECK(tagaspi_proc_init(GASPI_BLOCK));

    unsigned short g_rank, g_nranks;
    CHECK(gaspi_proc_rank(&g_rank));
    CHECK(gaspi_proc_num(&g_nranks));

    /* Should be the same as MPI */
    if (g_rank != sim->rank)
        die("wrong gaspi rank\n");

    if (g_nranks != sim->nranks)
        die("wrong gaspi nranks\n");

    /* Divide the segments into slots, one per box */

    /*
     * segment start                          segment end
     * |---------------------------------- ... ---------|
     * .                                                .
     * slot 0 (box 0)           slot 1 (box 1)  ...     .
     * |------------------------|--------- ...  --------|
     * .                        .                       .
     * .     pbbuf start         .                       .
     * |-----|------------------|-----|--- ... ---------|
     * | PBH |
     * |     |
     * |<--->|
     *  pboffset 
     */

    Gaspi *g = &sim->gaspi;
    g->nalloc = 16 * 1024; /* FIXME: compute */
    g->atomsize = NDIM * sizeof(double);
    g->pbsize = g->nalloc * g->atomsize;
    /* The slot contains the header at the beginning */
    g->slotsize = packbuf_data_size(g->nalloc, NDIM);
    g->pboffset = packbuf_data_size(0, 0);
    g->nslots = sim->nboxes;
    g->segsize = g->nslots * g->slotsize;

    for (int ineigh = 0; ineigh < NNEIGH; ineigh++) {
        for (int dir = 0; dir < PB_NDIR; dir++) {
            setup_gaspi_segment(sim, ineigh, dir, g->segsize);
        }
    }

    /* Setup queues */
    gaspi_number_t nqueues;
    CHECK(gaspi_queue_num(&nqueues));
    CHECK(tagaspi_queue_group_create(0, 0, nqueues,
                GASPI_QUEUE_GROUP_POLICY_CPU_RR));

    g->nqueues = nqueues;
}

static void
get_send_pair(Box *box, Neigh *neigh, enum pb_dir dir,
        int *sendibox, int *sendineigh)
{
    /* Always used the sending box/neigh as the reference value. In the
     * receiving part we compute the corresponding sending values */

    if (dir == PB_SEND) {
        *sendibox = box->i;
        *sendineigh = neigh->i;
    } else {
        *sendibox = neigh->boxid;
        *sendineigh = neigh->opposite->i;
    }
}

static int
build_gaspi_tag(int sendibox)
{
    return sendibox;
}

static void
setup_packbuf_gaspi(Sim *sim, Box *box, Neigh *neigh,
        PackBuf *pb, enum pb_type type, enum pb_dir dir)
{
    int sendibox, sendineigh;
    get_send_pair(box, neigh, dir, &sendibox, &sendineigh);

    /* Prevent accidental use */
    neigh = NULL;
    box = NULL;

    /* Use the send neigh index for both directions */
    int sendiseg = get_isegment(sendineigh, PB_SEND);
    int recviseg = get_isegment(sendineigh, PB_RECV);

    int tag = build_gaspi_tag(sendibox);

    size_t index = sendibox;
    size_t slot_offset = index * sim->gaspi.slotsize;
    int queue = sendibox % sim->gaspi.nqueues;

    /*
     * segment start                          segment end
     * .                                                .
     * seg                                              .
     * |---------------------------------- ... ---------|
     * .                                                .
     * slot 0 (box 0)           slot 1 (box 1)  ...     .
     * |------------------------|--------- ...  --------|
     * .                        .                       .
     * .                        .                       .
     * |-----|------------------|-----|--- ... ---------|
     *    ^                     |     |
     *    |                     |     pbbuf
     *    PackBufHeader         |
     *                          slot
     */

    void *seg = sim->gaspi.buf[sendineigh][dir];
    void *slot = seg + slot_offset;
    PackBufData *pbdata = slot;

    packbuf_gaspi_init(pb, pbdata,
            sendiseg, slot_offset,
            recviseg, slot_offset,
            sim->gaspi.nalloc, queue, tag);
}

static int
build_mpi_tag(int sendineigh, enum pb_req reqtype)
{
    int tag = sendineigh * PB_NREQS + reqtype;

    /* Ensure the tag is within the MPI standard limit */
    if (tag >= 32767) {
        die("tag exceed limit: %d >= %d\n", tag, 32767);
    }

    return tag;
}

static void
setup_packbuf_mpi(Sim *sim, Box *box, Neigh *neigh,
        PackBuf *pb, enum pb_type type, enum pb_dir dir)
{
    int sendibox, sendineigh;
    get_send_pair(box, neigh, dir, &sendibox, &sendineigh);

    int icomm = sendibox;
    MPI_Comm *comm = &sim->box[icomm].comm[type];

    int tags[PB_NREQS] = {
        [PB_BUF]    = build_mpi_tag(sendineigh, PB_BUF),
        [PB_NATOMS] = build_mpi_tag(sendineigh, PB_NATOMS)
    };

    packbuf_mpi_init(pb, tags, icomm, comm);
}

static void
setup_packbuf_comm(Sim *sim, Box *box, Neigh *neigh,
        PackBuf *pb, enum pb_type type, enum pb_dir dir)
{
    /* Enable GASPI only on the R packbuf */
    if (ENABLE_GASPI && type == PB_R) {
        setup_packbuf_gaspi(sim, box, neigh, pb, type, dir);
    } else {
        setup_packbuf_mpi(sim, box, neigh, pb, type, dir);
    }
}

static void
setup_packbuf_neigh(Sim *sim, Box *box, Neigh *neigh)
{
    /* Which PackBuf use the selection mechanism */
    int enablesel[PB_NTYPES] = {
        [PB_R]      = 0,
        [PB_RT]     = 1,
        [PB_RVT]    = 0
    };

    /* Setup the number of doubles needed per buffer */
    int ndoubles[PB_NTYPES] = {
        [PB_R]      = NDIM,
        [PB_RT]     = NDIM + 1,
        [PB_RVT]    = NDIM + NDIM + 1
    };

    int remoterank[PB_NDIR] = {
        [PB_SEND] = neigh->rank,
        [PB_RECV] = neigh->rank
    };

    for (int type = 0; type < PB_NTYPES; type++) {
        for (int dir = 0; dir < PB_NDIR; dir++) {
            PackBuf *pb = &neigh->pb[type][dir];
            packbuf_init(pb, enablesel[type], ndoubles[type], remoterank[dir]);

            sprintf(pb->name, "PackBuf[type=%s dir=%s rank=%d box=%d neigh=%d]",
                    PB_TYPENAME(type), PB_DIRNAME(dir),
                    sim->rank, box->i, neigh->i);

            setup_packbuf_comm(sim, box, neigh, pb, type, dir);
            packbuf_debug_switch(pb, PB_GARBAGE, PB_READY);

            /* Set the PB pointers in the box table */
            box->pb[neigh->i][type][dir] = pb;
        }
    }
}

static void
setup_mpi(Sim *sim)
{
    /* Use one communicator per box */
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < PB_NTYPES; j++) {
            MPI_Comm_dup(MPI_COMM_WORLD, &box->comm[j]);
        }
    }
}

void
setup_packbuf(Sim *sim)
{
    setup_mpi(sim);

    if (ENABLE_GASPI)
        setup_gaspi_segments(sim);

    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            setup_packbuf_neigh(sim, box, &box->neigh[j]);
        }
    }
}
