#ifndef COMM_H
#define COMM_H

#include "types.h"

void comm_pack  (Sim *sim, enum pb_type type, enum pb_dir dir);
void comm_unpack(Sim *sim, enum pb_type type, enum pb_dir dir);

void comm_send(Sim *sim, enum pb_type type, enum pb_dir dir, enum pb_req req);
void comm_recv(Sim *sim, enum pb_type type, enum pb_dir dir, enum pb_req req);
void comm_wait(Sim *sim, enum pb_type type, enum pb_dir dir, enum pb_req req);

void comm_signal(Sim *sim, enum pb_type type, enum pb_dir dir);
void comm_linger(Sim *sim, enum pb_type type, enum pb_dir dir);

void comm_setup(Sim *sim);
void comm_tidy(Sim *sim);
void comm_borders(Sim *sim);
void comm_ghost_position(Sim *sim);
void comm_waitall(Sim *sim);

void comm_ready(Sim *sim);
void comm_check_header(Sim *sim);
void comm_reset_header(Sim *sim);

#endif /* COMM_H */
