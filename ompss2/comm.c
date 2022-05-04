#include "types.h"
#include "comm.h"

/** Transfers local atoms outside the box domain to the correct box */
void
comm_tidy(Sim *sim)
{
    enum pb_type type = PB_RVT;

    comm_wait  (sim, type, PB_SEND, PB_NATOMS);
    comm_wait  (sim, type, PB_SEND, PB_BUF);
    comm_pack  (sim, type, PB_SEND);
    comm_send  (sim, type, PB_SEND, PB_NATOMS);
    comm_send  (sim, type, PB_SEND, PB_BUF);

    comm_recv  (sim, type, PB_RECV, PB_NATOMS);
    comm_wait  (sim, type, PB_RECV, PB_NATOMS);
    comm_recv  (sim, type, PB_RECV, PB_BUF);
    comm_wait  (sim, type, PB_RECV, PB_BUF);
    comm_unpack(sim, type, PB_RECV);
}

/** Transfers the information (count, position and type) of local atoms in the
 * borders to the ghost atoms in the neighboring boxes. */
void
comm_borders(Sim *sim)
{
    enum pb_type type = PB_RT;

    comm_wait  (sim, type, PB_SEND, PB_NATOMS);
    comm_wait  (sim, type, PB_SEND, PB_BUF);
    comm_pack  (sim, type, PB_SEND);
    comm_send  (sim, type, PB_SEND, PB_NATOMS);
    comm_send  (sim, type, PB_SEND, PB_BUF);

    comm_recv  (sim, type, PB_RECV, PB_NATOMS);
    comm_wait  (sim, type, PB_RECV, PB_NATOMS);
    comm_recv  (sim, type, PB_RECV, PB_BUF);
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

    comm_wait  (sim, type, PB_SEND, PB_BUF);
    comm_pack  (sim, type, PB_SEND);
    comm_send  (sim, type, PB_SEND, PB_BUF);

    comm_recv  (sim, type, PB_RECV, PB_BUF);
    comm_wait  (sim, type, PB_RECV, PB_BUF);
    comm_unpack(sim, type, PB_RECV);
}

/** Waits for all communications to finish. */
void
comm_waitall(Sim *sim)
{
    for (int i = 0; i < PB_NTYPES; i++) {
        for (int j = 0; j < PB_NDIR; j++) {
            comm_wait(sim, i, j, PB_BUF);
            comm_wait(sim, i, j, PB_NATOMS);
        }
    }

    #pragma oss taskwait /* required */
}
