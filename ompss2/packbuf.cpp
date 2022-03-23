//#define ENABLE_DEBUG
#include "log.h"
#include "types.h"

#include <stdlib.h>
#include <mpi.h>

#define PACKBUF_INCR 2000

static void
packbuf_switch(PackBuf *pb, enum packbuf_state prev, enum packbuf_state next)
{
    if (pb->state != prev)
        die("packbuf in state %d, expected %d\n", pb->state, prev);

    pb->state = next;
}

void
packbuf_debug_switch(PackBuf *pb, enum packbuf_state prev, enum packbuf_state next)
{
    if (pb->debug_state != prev)
        die("packbuf in debug_state %d, expected %d\n", pb->debug_state, prev);

    pb->debug_state = next;
}

/* Grows the buffer so that the allocated capacity can hold at least n
 * atoms. */
static void
packbuf_grow(PackBuf *pb, int n)
{
    if (n < pb->natoms)
        n = pb->natoms;

    if (pb->nalloc < n) {
        //if (pb->nalloc + PACKBUF_INCR < n)
        //    n = pb->nalloc + PACKBUF_INCR;

        size_t newalloc = sizeof(double) * pb->atomsize * n;

        if (newalloc == 0)
            abort();

        pb->buf = (double *) safe_realloc(pb->buf, newalloc);

        /* Also grow the selection buffer if enabled */
        if (pb->enable_sel) {
            pb->sel = (int *) safe_realloc(pb->sel, sizeof(int) * n);
        }

        pb->nalloc = n;
    }
}

/* Ensures the buffer can hold at least nextra additional atoms */
static void
packbuf_grow_extra(PackBuf *pb, int nextra)
{
    packbuf_grow(pb, pb->natoms + nextra);
}

void
packbuf_mpisend_buf(PackBuf *pb)
{
    packbuf_switch(pb, PB_READY, PB_SENDING);

    dbg("packbuf_mpisend_buf: natoms=%d remoterank=%d tag=%d\n",
            pb->natoms, pb->remoterank, pb->tag);

	if (ENABLE_NONBLOCKING_MPI) {
		if (pb->waitreq)
			die("packbuf_mpisend_buf: buffer in use\n");

		if (pb->natoms != 0) {
			MPI_Isend((void *) pb->buf, pb->natoms * pb->atomsize,
					MPI_DOUBLE, pb->remoterank, pb->tag, *pb->comm, &pb->req);
			pb->waitreq = 1;
		}
	} else {
		if (pb->natoms != 0) {
			MPI_Send((void *) pb->buf, pb->natoms * pb->atomsize,
					MPI_DOUBLE, pb->remoterank, pb->tag, *pb->comm);
		}
	}

    packbuf_switch(pb, PB_SENDING, PB_READY);
}

void
packbuf_mpisend(PackBuf *pb)
{
    dbg("packbuf_mpisend: natoms=%d remoterank=%d tag=%d\n",
            pb->natoms, pb->remoterank, pb->tag);

	if (ENABLE_NONBLOCKING_MPI && pb->waitreqn)
		die("packbuf_mpisend: buffer in use\n");

	void *buf = (void *) &pb->natoms;

	if (ENABLE_NONBLOCKING_MPI) {
		MPI_Isend(buf, 1, MPI_INT, pb->remoterank, pb->tag, *pb->comm,
				&pb->reqn);

		pb->waitreqn = 1;
	} else {
		MPI_Send(buf, 1, MPI_INT, pb->remoterank, pb->tag, *pb->comm);
	}

    packbuf_mpisend_buf(pb);
}

void
packbuf_mpirecv_buf(PackBuf *pb, int natoms)
{
    packbuf_switch(pb, PB_READY, PB_RECVING);

    dbg("packbuf_mpirecv_buf: natoms=%d remoterank=%d tag=%d\n",
            natoms, pb->remoterank, pb->tag);

	if (ENABLE_NONBLOCKING_MPI && pb->waitreq)
		die("packbuf_mpirecv_buf: buffer in use\n");

	if (natoms > 0) {
		/* Grow the buffer if needed */
		packbuf_grow(pb, natoms);

		/* And receive that many atoms */
		int size = natoms * pb->atomsize;

		if (ENABLE_NONBLOCKING_MPI) {
			MPI_Irecv((void *) pb->buf, size, MPI_DOUBLE,
					pb->remoterank, pb->tag, *pb->comm, &pb->req);
			pb->waitreq = 1;
		} else {
			MPI_Recv((void *) pb->buf, size, MPI_DOUBLE,
					pb->remoterank, pb->tag, *pb->comm, MPI_STATUS_IGNORE);
		}
	}

	/* FIXME: this is dangerous as we are writing the natoms in the buffer
	while they may be still being written by MPI_Irecv */
    pb->natoms = natoms;
    packbuf_switch(pb, PB_RECVING, PB_READY);
}

void
packbuf_mpirecv_natoms(PackBuf *pb)
{
    /* Find out how many atoms I need to make room for */
    dbg("packbuf_mpirecv_natoms: natoms=? remoterank=%d tag=%d\n",
            pb->remoterank, pb->tag);

	if (ENABLE_NONBLOCKING_MPI) {
		if (pb->waitreqn)
			die("packbuf_mpirecv_natoms: buffer in use\n");

		MPI_Irecv((void *) &pb->recvnatoms, 1, MPI_INT,
				pb->remoterank, pb->tag, *pb->comm, &pb->reqn);

		pb->waitreqn = 1;
	} else {
		MPI_Recv((void *) &pb->recvnatoms, 1, MPI_INT,
				pb->remoterank, pb->tag, *pb->comm, MPI_STATUS_IGNORE);
	}
}

void
packbuf_shmcopy(PackBuf *src, PackBuf *dst)
{
    packbuf_switch(src, PB_READY, PB_COPYING);
    packbuf_switch(dst, PB_READY, PB_COPYING);

    if (src->natoms != 0) {
        packbuf_grow(dst, src->natoms);
        memcpy(dst->buf, src->buf, src->natoms * src->atomsize * sizeof(double));
    }
    dst->natoms = src->natoms;

    packbuf_switch(dst, PB_COPYING, PB_READY);
    packbuf_switch(src, PB_COPYING, PB_READY);
}

/* FIXME: move to .h so the compiler can optimize the constant NULL
 * pointers */
void
packbuf_add(PackBuf *pb, Vec *r, Vec *v, int *type)
{
    packbuf_switch(pb, PB_READY, PB_ADDING);
    /* Ensure we have room for another atom */
    packbuf_grow_extra(pb, 1);

    if (pb->waitreqn || pb->waitreq)
        die("packbuf_add: buffer in use\n");

    int j = pb->natoms * pb->atomsize;

    if (r != NULL) {
        for (int d = X; d <= Z; d++)
            pb->buf[j++] = r[0][d];
    }

    if (v != NULL) {
        for (int d = X; d <= Z; d++)
            pb->buf[j++] = v[0][d];
    }

    if (type != NULL) {
        /* FIXME: We are sending the type as a double */
        pb->buf[j++] = (double) type[0];
    }

    pb->natoms++;

//    dbg("packbuf %p now has %d atoms (alloc %d)\n",
//            pb, pb->natoms, pb->nalloc);

    if (j != pb->natoms * pb->atomsize)
        die("packbuf_add atom size mismatch\n");

    packbuf_switch(pb, PB_ADDING, PB_READY);
}

void
packbuf_add_sel(PackBuf *pb, Vec *r, Vec *v, int *type, int iatom)
{
    int ncur = pb->natoms;
    packbuf_add(pb, r, v, type);
    pb->sel[ncur] = iatom;
}

void
packbuf_unpack(PackBuf *pb, Vec *r, Vec *v, int *types)
{
    packbuf_switch(pb, PB_READY, PB_UNPACKING);

    if (pb->waitreqn || pb->waitreq)
        die("packbuf_unpack: buffer in use\n");

    for (int i = 0, j = 0; i < pb->natoms; i++) {
        if (r != NULL) {
            for (int d = X; d <= Z; d++)
                r[i][d] = pb->buf[j++];
        }

        if (v != NULL) {
            for (int d = X; d <= Z; d++)
                v[i][d] = pb->buf[j++];
        }

        if (types != NULL) {
            types[i] = (int) pb->buf[j++];
        }
    }
    packbuf_switch(pb, PB_UNPACKING, PB_READY);
}

void
packbuf_unpack_sel(PackBuf *pb, Vec *r, Vec *v, int *types, int *sel)
{
    packbuf_switch(pb, PB_READY, PB_UNPACKING);

    if (pb->waitreqn || pb->waitreq)
        die("packbuf_unpack_sel: buffer in use\n");

    for (int i = 0, j = 0; i < pb->natoms; i++) {
        if (r != NULL) {
            for (int d = X; d <= Z; d++)
                r[sel[i]][d] = pb->buf[j++];
        }

        if (v != NULL) {
            for (int d = X; d <= Z; d++)
                v[sel[i]][d] = pb->buf[j++];
        }

        if (types != NULL) {
            types[sel[i]] = (int) pb->buf[j++];
        }
    }
    packbuf_switch(pb, PB_UNPACKING, PB_READY);
}

void
packbuf_clear(PackBuf *pb)
{
    packbuf_switch(pb, PB_READY, PB_CLEANING);

    if (pb->waitreqn || pb->waitreq)
        die("packbuf_unpack_sel: buffer in use\n");

    pb->natoms = 0;

    packbuf_switch(pb, PB_CLEANING, PB_READY);
}

void
packbuf_init(PackBuf *pb, int enable_sel, int atomsize,
        int remoterank, int tag, MPI_Comm *comm)
{
    memset(pb, 0, sizeof(*pb));

    pb->atomsize = atomsize;
    pb->enable_sel = enable_sel;
    pb->comm = comm;

    if (remoterank < 0)
        abort();

    pb->remoterank = remoterank;
    pb->tag = tag;
    packbuf_switch(pb, PB_GARBAGE, PB_READY);
}
