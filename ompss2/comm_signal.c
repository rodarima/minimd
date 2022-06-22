#define ENABLE_DEBUG 0
#include "comm.h"
#include "types.h"
#include "log.h"
#include "packbuf.h"

void
comm_signal(Sim *sim, enum pb_type type, enum pb_dir dir)
{
    for (int ibox = 0; ibox < sim->nboxes; ibox++) {
        Box *box = &sim->box[ibox];
        for (int ineigh = 0; ineigh < NNEIGH; ineigh++) {
            Neigh *neigh = &box->neigh[ineigh];
            PackBuf *pb = &neigh->pb[type][dir];

            #pragma oss task label("comm_signal") \
                inout(pb->natoms) inout(pb->data) firstprivate(pb)
            packbuf_signal(pb);
        }
    }
}

void
comm_linger(Sim *sim, enum pb_type type, enum pb_dir dir)
{
    for (int ibox = 0; ibox < sim->nboxes; ibox++) {
        Box *box = &sim->box[ibox];
        for (int ineigh = 0; ineigh < NNEIGH; ineigh++) {
            Neigh *neigh = &box->neigh[ineigh];
            PackBuf *pb = &neigh->pb[type][dir];

            #pragma oss task label("comm_linger") \
                inout(pb->natoms) inout(pb->data) firstprivate(pb)
            packbuf_linger(pb);
        }
    }
}

void
comm_check_header(Sim *sim)
{
    for (int ibox = 0; ibox < sim->nboxes; ibox++) {
        Box *box = &sim->box[ibox];
        for (int ineigh = 0; ineigh < NNEIGH; ineigh++) {
            Neigh *neigh = &box->neigh[ineigh];
            for (int type = 0; type < PB_NTYPES; type++) {
                for (int dir = 0; dir < PB_NDIR; dir++) {
                    PackBuf *pb = &neigh->pb[type][dir];

                    packbuf_header_check(pb);
                }
            }
        }
    }
}

void
comm_reset_header(Sim *sim)
{
    for (int ibox = 0; ibox < sim->nboxes; ibox++) {
        Box *box = &sim->box[ibox];
        for (int ineigh = 0; ineigh < NNEIGH; ineigh++) {
            Neigh *neigh = &box->neigh[ineigh];
            for (int type = 0; type < PB_NTYPES; type++) {
                for (int dir = 0; dir < PB_NDIR; dir++) {
                    PackBuf *pb = &neigh->pb[type][dir];

                    packbuf_header_reset(pb);
                }
            }
        }
    }
}
