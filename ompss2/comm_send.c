#define ENABLE_DEBUG 1
#include "types.h"
#include "log.h"
#include "packbuf.h"

#pragma oss task label("box_send_neigh_task") \
    in(pb->buf) in(pb->natoms)
static void
box_send_neigh_task(Sim *sim, Box *box, Neigh *neigh, PackBuf *pb,
        enum pb_type type)
{
    packbuf_debug_switch(pb, PB_READY, PB_SENDING);

    if (type == PB_SEND_R) {
        if (ENABLE_GASPI) {
            packbuf_gaspi_send_buf(pb);
        } else {
            packbuf_mpi_send_buf(pb);
        }

    } else {
        packbuf_mpi_send(pb);
    }

    packbuf_debug_switch(pb, PB_SENDING, PB_READY);
}

static void
box_send_neigh(Sim *sim, Box *box, Neigh *neigh, enum pb_type type)
{
    if (neigh->rank == sim->rank) {
        /* No-op: will be copied at recv */
        return;
    }

    PackBuf *pb = box->pb[type][neigh->i];

    box_send_neigh_task(sim, box, neigh, pb, type);
}

void
comm_send(Sim *sim, enum pb_type type)
{
    dbg("comm_send %s\n", PB_TYPENAME(type));
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            box_send_neigh(sim, box, &box->neigh[j], type);
        }
    }

    #pragma oss taskwait
}
