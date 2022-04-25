#ifndef PACKBUF_H
#define PACKBUF_H

#include <mpi.h>
#include <stdlib.h>

void packbuf_init(PackBuf *pb,
        int enable_sel, int atomsize, int remoterank,
        int tag[PB_NREQTYPES], int icomm, MPI_Comm *comm);

void packbuf_switch(PackBuf *pb, enum packbuf_state prev, enum packbuf_state next);
void packbuf_debug_switch(PackBuf *pb, enum packbuf_state prev, enum packbuf_state next);
void packbuf_shmcopy(PackBuf *src, PackBuf *dst, enum pb_reqtype reqtype);
void packbuf_add(PackBuf *pb, Vec *r, Vec *v, int *type);
void packbuf_add_sel(PackBuf *pb, Vec *r, Vec *v, int *type, int iatom);
void packbuf_unpack(PackBuf *pb, Vec *r, Vec *v, int *types);
void packbuf_unpack_sel(PackBuf *pb, Vec *r, Vec *v, int *types, int *sel);
void packbuf_clear(PackBuf *pb);
void packbuf_grow(PackBuf *pb, int n);

void packbuf_send(PackBuf *pb, enum pb_reqtype reqtype);
void packbuf_recv(PackBuf *pb, enum pb_reqtype reqtype);

enum pb_type packbuf_opposite_dir(enum pb_type type);

/* MPI */

void packbuf_mpi_send(PackBuf *pb, enum pb_reqtype reqtype);
void packbuf_mpi_recv(PackBuf *pb, enum pb_reqtype reqtype);
void packbuf_mpi_waitn(PackBuf **pbs, int n, enum pb_reqtype reqtype);

/* GASPI */

void packbuf_gaspi_init(PackBuf *pb, double *newbuf,
        int sendseg, size_t send_offset_bytes,
        int recvseg, size_t recv_offset_bytes,
        size_t nalloc, int queue);

void packbuf_gaspi_send(PackBuf *pb, enum pb_reqtype reqtype);
void packbuf_gaspi_recv(PackBuf *pb, enum pb_reqtype reqtype);

#endif /* PACKBUF_H */
