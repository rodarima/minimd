#include "types.h"
#include "comm.h"

void
comm_tidy(Sim *sim)
{
    comm_wait(sim, PB_SEND_RVT, PB_NATOMS);
    comm_wait(sim, PB_SEND_RVT, PB_BUF);
    comm_pack(sim, PB_SEND_RVT);
    comm_send(sim, PB_SEND_RVT);

    comm_recv(sim, PB_RECV_RVT, PB_NATOMS);
    comm_wait(sim, PB_RECV_RVT, PB_NATOMS);
    comm_recv(sim, PB_RECV_RVT, PB_BUF);
    comm_wait(sim, PB_RECV_RVT, PB_BUF);
    comm_unpack(sim, PB_RECV_RVT);
}

void
comm_borders(Sim *sim)
{
    comm_wait(sim, PB_SEND_RT, PB_BUF);
    comm_wait(sim, PB_SEND_RT, PB_NATOMS);
    comm_pack(sim, PB_SEND_RT);
    comm_send(sim, PB_SEND_RT);

    comm_recv(sim, PB_RECV_RT, PB_NATOMS);
    comm_wait(sim, PB_RECV_RT, PB_NATOMS);
    comm_recv(sim, PB_RECV_RT, PB_BUF);
    comm_wait(sim, PB_RECV_RT, PB_BUF);
    comm_unpack(sim, PB_RECV_RT);
}

void
comm_ghost_position(Sim *sim)
{
    comm_wait(sim, PB_SEND_R, PB_BUF);
    comm_pack(sim, PB_SEND_R);
#pragma oss taskwait
    comm_send(sim, PB_SEND_R);

    comm_recv(sim, PB_RECV_R, PB_BUF);
    comm_wait(sim, PB_RECV_R, PB_BUF);
    comm_unpack(sim, PB_RECV_R);
}
