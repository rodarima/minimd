#ifndef COMM_H
#define COMM_H

#include "types.h"

void comm_pack  (Sim *sim, enum pb_type type, enum pb_dir dir);
void comm_unpack(Sim *sim, enum pb_type type, enum pb_dir dir);

void comm_send(Sim *sim, enum pb_type type, enum pb_dir dir, enum pb_req req);
void comm_recv(Sim *sim, enum pb_type type, enum pb_dir dir, enum pb_req req);
void comm_wait(Sim *sim, enum pb_type type, enum pb_dir dir, enum pb_req req);

void comm_setup(Sim *sim);
void comm_tidy(Sim *sim);
void comm_borders(Sim *sim);
void comm_ghost_position(Sim *sim);
void comm_waitall(Sim *sim);


#endif /* COMM_H */
