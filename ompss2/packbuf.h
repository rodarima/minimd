#ifndef PACKBUF_H
#define PACKBUF_H

#include <mpi.h>
#include <stdlib.h>

/* Used to mark the status of the buffer, so we can detect concurrent
 * access to the buffer and explain what was using it before */
enum packbuf_state {
    PB_GARBAGE = 0,
    PB_READY = 1,
    PB_PACKING = 2,
    PB_SENDING = 3,
    PB_RECVING = 4,
    PB_COPYING = 5,
    PB_UNPACKING = 6,
    PB_ADDING = 7,
    PB_READING = 8,
    PB_CLEANING = 9,
    PB_WAITING = 10,
};

enum pb_magic {
    PB_MAGIC_OK = 12345,
    PB_MAGIC_KO = 666
};

enum gaspi_segment_dir {
    SENDSEG = 0,
    RECVSEG = 1
};

enum pb_req {
    PB_BUF = 0,
    PB_NATOMS = 1,
    PB_NREQS = 2
};

enum pb_mode {
    PB_BAD = 0,
    PB_MPI = 1,
    PB_GASPI = 2
};

enum pb_type {
    PB_R        = 0,
    PB_RT       = 1,
    PB_RVT      = 2,
    PB_NTYPES   = 3
};

enum pb_dir {
    PB_SEND = 0,
    PB_RECV = 1,
    PB_NDIR = 2
};

/* Transform enum constants into names */
#define PB_TYPENAME(x)  (((char *[]){"R","RT","RVT"})[x])
#define PB_DIRNAME(x)   (((char *[]){"SEND","RECV"})[x])
#define PB_REQNAME(x)   (((char *[]){"BUF","NATOMS"})[x])

typedef struct {
    int tag[PB_NREQS];
    MPI_Request req[PB_NREQS];
    MPI_Comm *comm;
    int icomm;                  /* And index to identify the MPI_Comm */
} PackBufMPI;

typedef struct {
    int nid; /* Notification id */
    int sendseg;
    int recvseg;
    size_t sendoffset; /* in bytes */
    size_t recvoffset;
    int queue;
} PackBufGASPI;

typedef struct {
    int magic;
    int srcbox;
    int senddir;
    int dstbox;
    int icomm;
} PackBufHeader;

/* The data to be communicated */
typedef struct {
    PackBufHeader header;
    int xnatoms;    /* Exchanged number of atoms */
    double buf[];   /* Exchanged floating point data */
} PackBufData;

typedef struct {
    int box;
    int senddir;

    enum pb_dir dir;

    char name[256];
    int natoms;     /* Number of atoms currently in the buffer */
    int nalloc;     /* Number of atoms allocated */
    int atomsize;   /* Number of doubles required per atom */

    size_t datasize;
    PackBufData *data;

    int enable_sel; /* If non-zero use selection for packing */
    int *sel;       /* Selection of atoms */

    int remoterank;
    int waitreq[PB_NREQS];  /* Wait before overwrite buf/natoms */

    enum packbuf_state debug_state; /* Reserved for debugging purposes */
    enum packbuf_state state;

    enum pb_mode mode;
    union {
        PackBufMPI mpi;
        PackBufGASPI gaspi;
    };
} PackBuf;

void packbuf_init(PackBuf *pb, enum pb_dir dir,
        int enable_sel, int atomsize, int remoterank);

void packbuf_switch(PackBuf *pb, enum packbuf_state prev, enum packbuf_state next);
void packbuf_debug_switch(PackBuf *pb, enum packbuf_state prev, enum packbuf_state next);
void packbuf_shmcopy(PackBuf *src, PackBuf *dst, enum pb_req req);
void packbuf_add(PackBuf *pb, Vec *r, Vec *v, int *type);
void packbuf_add_sel(PackBuf *pb, Vec *r, Vec *v, int *type, int iatom);
void packbuf_unpack(PackBuf *pb, Vec *r, Vec *v, int *types);
void packbuf_unpack_sel(PackBuf *pb, Vec *r, Vec *v, int *types, int *sel);
void packbuf_clear(PackBuf *pb);
void packbuf_grow(PackBuf *pb, int n);

void packbuf_send(PackBuf *pb, enum pb_req req);
void packbuf_recv(PackBuf *pb, enum pb_req req);

size_t packbuf_data_size(size_t natoms, size_t atomdoubles);

void packbuf_check_header(PackBuf *pb);

/* MPI */

void packbuf_mpi_init(PackBuf *pb, int tag[PB_NREQS],
        int icomm, MPI_Comm *comm);

void packbuf_mpi_send(PackBuf *pb, enum pb_req req);
void packbuf_mpi_recv(PackBuf *pb, enum pb_req req);
void packbuf_mpi_waitn(PackBuf **pbs, int n, enum pb_req req);

/* GASPI */

void packbuf_gaspi_init(PackBuf *pb, PackBufData *newdata,
        int sendseg, size_t send_offset_bytes,
        int recvseg, size_t recv_offset_bytes,
        size_t nalloc, int queue, int tag);

void packbuf_gaspi_send(PackBuf *pb, enum pb_req req);
void packbuf_gaspi_recv(PackBuf *pb, enum pb_req req);

#endif /* PACKBUF_H */
