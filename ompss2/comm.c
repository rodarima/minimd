#define ENABLE_DEBUG 0
#include "log.h"
#include "types.h"
#include "comm.h"
#include "trace.h"

/** Transfers local atoms outside the box domain to the correct box */
void
comm_tidy(Sim *sim)
{
    enum pb_type type = PB_RVT;

    comm_pack  (sim, type, PB_SEND);

    comm_send  (sim, type, PB_SEND, PB_NATOMS);
    comm_recv  (sim, type, PB_RECV, PB_NATOMS);
    comm_wait  (sim, type, PB_SEND, PB_NATOMS);
    comm_wait  (sim, type, PB_RECV, PB_NATOMS);

    comm_send  (sim, type, PB_SEND, PB_BUF);
    comm_recv  (sim, type, PB_RECV, PB_BUF);
    comm_wait  (sim, type, PB_SEND, PB_BUF);
    comm_wait  (sim, type, PB_RECV, PB_BUF);

    comm_unpack(sim, type, PB_RECV);
}

static void
fence(Sim *sim, const char *name)
{
    #pragma oss taskwait
    MPI_Barrier(MPI_COMM_WORLD);
    sleep(1);
    MPI_Barrier(MPI_COMM_WORLD);
    err("rank=%d: COMM BORDERS: %s OK\n", sim->rank, name);
    MPI_Barrier(MPI_COMM_WORLD);
    sleep(1);
    MPI_Barrier(MPI_COMM_WORLD);

}

/** Transfers the information (count, position and type) of local atoms in the
 * borders to the ghost atoms in the neighboring boxes. */
void
comm_borders(Sim *sim)
{
    enum pb_type type = PB_RT;

    comm_pack  (sim, type, PB_SEND);

    comm_send  (sim, type, PB_SEND, PB_NATOMS);
    comm_recv  (sim, type, PB_RECV, PB_NATOMS);

    comm_wait  (sim, type, PB_SEND, PB_NATOMS);
    comm_wait  (sim, type, PB_RECV, PB_NATOMS);

    comm_send  (sim, type, PB_SEND, PB_BUF);
    comm_recv  (sim, type, PB_RECV, PB_BUF);

    comm_wait  (sim, type, PB_SEND, PB_BUF);
    comm_wait  (sim, type, PB_RECV, PB_BUF);

    comm_unpack(sim, type, PB_RECV);
}

/** Transfers the updated positions of the local atoms in the borders. The
 * number of atoms is always constant between iterations, even if they have
 * moved to another box. */
void
comm_ghost_position(Sim *sim)
{
    enum pb_type type = PB_R;

    comm_pack  (sim, type, PB_SEND);
    comm_linger(sim, type, PB_SEND);
    comm_send  (sim, type, PB_SEND, PB_BUF);
    comm_recv  (sim, type, PB_RECV, PB_BUF);
    comm_wait  (sim, type, PB_SEND, PB_BUF);
    comm_wait  (sim, type, PB_RECV, PB_BUF);
    comm_unpack(sim, type, PB_RECV);
    comm_signal(sim, type, PB_RECV);
}

/** Waits for all communications to finish. */
void
comm_waitall(Sim *sim)
{
    for (int i = 0; i < PB_NTYPES; i++) {
        for (int j = 0; j < PB_NDIR; j++) {
            comm_wait(sim, i, j, PB_BUF);
            comm_wait(sim, i, j, PB_NATOMS);
            #pragma oss taskwait /* required */
        }
    }
}

void
comm_ready(Sim *sim)
{
    comm_signal(sim, PB_R, PB_RECV);
}
