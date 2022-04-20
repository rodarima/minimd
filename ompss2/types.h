/* ----------------------------------------------------------------------
   miniMD is a simple, parallel molecular dynamics (MD) code.   miniMD is
   an MD microapplication in the Mantevo project at Sandia National
   Laboratories ( http://www.mantevo.org ). The primary
   authors of miniMD are Steve Plimpton (sjplimp@sandia.gov) , Paul Crozier
   (pscrozi@sandia.gov) and Christian Trott (crtrott@sandia.gov).

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This library is free software; you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation;
   either version 3 of the License, or (at your option) any later
   version.

   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this software; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA.  See also: http://www.gnu.org/licenses/lgpl.txt .

   For questions, contact Paul S. Crozier (pscrozi@sandia.gov) or
   Christian Trott (crtrott@sandia.gov).

   Please read the accompanying README and LICENSE files.
---------------------------------------------------------------------- */

#ifndef TYPES_H
#define TYPES_H

#include <mpi.h>
#include <stdlib.h>

#define VARIANT_STRING "miniMD-ompss 2.0 (MPI+OmpSs-2+TAMPI)"

enum ForceStyle { FORCELJ, FORCEEAM };

#ifndef PAD4
#define PAD 3
#else
#define PAD 4
#endif

#ifdef __INTEL_COMPILER
#ifndef ALIGNMALLOC
#define ALIGNMALLOC 64
#endif
#define RESTRICT __restrict
#endif

#ifndef RESTRICT
#define RESTRICT
#endif

/* Dimensions */
#define X 0
#define Y 1
#define Z 2
#define NDIM 3

/* Limits per dimension */
#define LO 0
#define HI 1
#define NLIM 2

typedef double Vec[PAD];
typedef double VecPack[NDIM];
typedef double Domain[NDIM][NLIM];
typedef int    Range[NDIM][NLIM];

#include "config.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

typedef struct sim Sim;

typedef struct bin {
    int natoms; /* Number of atoms present in this bin */
    int nalloc; /* Number of atoms allocated in the bin */
    int *atom; /* Indexes of the atoms */

    /* Thermal information keep per time step in reduced units */
    double pot_energy; /* Potential energy */
    double kin_energy; /* Kinetic energy */
    double vdwl_energy; /* Van der Waals pairwise energy */
    double virial_pressure; /* Virial temperature */
    double potghost_energy; /* Potential energy of the ghosts only */
} Bin;

#define NNEIGHSIDE 1
#define NNEIGHDIM (NNEIGHSIDE*2 + 1) 
#define NNEIGH (NNEIGHDIM*NNEIGHDIM*NNEIGHDIM - 1)
#define NSUB 27

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

enum gaspi_segment_dir {
    SENDSEG = 0,
    RECVSEG = 1
};

enum pb_reqtype {
    PB_BUF = 0,
    PB_NATOMS = 1,
    PB_NREQTYPES = 2
};

typedef struct {
    int natoms;     /* Number of atoms currently in the buffer */
    int recvnatoms; /* Number of atoms to be received */
    int nalloc;     /* Number of atoms allocated */
    int atomsize;   /* Number of doubles required per atom */
    double *buf;    /* The contiguous buffer */
    int enable_sel; /* If non-zero use selection for packing */
    int *sel;       /* Selection of atoms */
    enum packbuf_state debug_state; /* Reserved for debugging purposes */
    enum packbuf_state state;
    MPI_Request req[PB_NREQTYPES];
    MPI_Comm *comm;
    int icomm;      /* And index to identify the MPI_Comm */
    int waitreq[PB_NREQTYPES];    /* Wait for the request before writing the buffer */
    int remoterank;
    int tag;

    /* GASPI related */
    int gaspi;
    int sendseg;
    int recvseg;
    size_t sendoffset; /* in bytes */
    size_t recvoffset;
    int queue;
} PackBuf;

typedef struct box Box;
typedef struct neigh Neigh;

/* A neighboring box */
typedef struct neigh {
    int i;              /* Local index of this neighbor in the box */
    int boxcoord[NDIM]; /* Corresponding box coordinate without wrapping */
    int boxcoordw[NDIM];/* Corresponding box coordinate wrapped */
    int boxid;          /* Corresponding box id */
    int delta[NDIM];    /* Delta vector in boxes */
    int rank;           /* Neighbor process rank */
    int rankcoord[NDIM];/* Neighbor process coordinates without wrapping */

    Box *box;           /* Neighbor box or NULL if outside the rank */
    Neigh *opposite;    /* Opposite neighbor at -delta */

    int wraps;  /* Non-zero if some dimension needs PBC correction */
    Vec addpbc; /* PBC correction per dimension */

    /* Communication packing buffers */
    PackBuf send_r;     /* Send atom positions */
    PackBuf recv_r;     /* Receive atom positions */
    PackBuf send_rt;    /* Send atom positions and types */
    PackBuf recv_rt;    /* Receive atom positions and types */
    PackBuf send_rvt;   /* Send atom position, velocity and type */
    PackBuf recv_rvt;   /* Receive atom positions, velocity and type */

} Neigh;

/* Subdivision of the box into cubic subdomains, each with a list of
 * neighboring boxes */
typedef struct Subdomain {
    int i;                  /* Subdomain index */
    int nneigh;             /* # of neighbors of this subdomain */
    int delta[NDIM];        /* Delta offset of the subdomain */
    Neigh *neigh[NNEIGH];   /* List of neighbors (pointers) */
} Subdomain;

/* Holds the list of nearby atoms (also known as neighbor atoms) */
typedef struct Nearby {
    int natoms;
    int nalloc;
    int *atom;
} Nearby;

#define HIST_MAX_NBINS 400

/* Histogram structure */
typedef struct Hist {
    int nbins;
    int count[HIST_MAX_NBINS];
    double delta;
    const char *filepath;
} Hist;

enum pb_type {
    PB_SEND_R = 0,
    PB_RECV_R,
    PB_SEND_RT,
    PB_RECV_RT,
    PB_SEND_RVT,
    PB_RECV_RVT,
    PB_NTYPES
};

#define PB_TYPENAME(x) ( \
((char *[]){ \
"PB_SEND_R", "PB_RECV_R", \
"PB_SEND_RT", "PB_RECV_RT", \
"PB_SEND_RVT", "PB_RECV_RVT" \
})[x])

#define PB_REQTYPENAME(x) ( \
((char *[]){ \
"PB_BUF", "PB_NATOMS" \
})[x])


/* All information needed for a box of the simulation */
typedef struct box {
    int i;          /* Box index for this process */
    int idim[NDIM]; /* Box index per dimension */

    int nlocal; /* Number of local atoms in this box */
    int nghost; /* Number of ghost atoms in this box */
    int nalloc; /* Allocated atom capacity including ghosts */

    int nbinscore[NDIM]; /* # of bins in the core domain */
    int nbinshalo[NDIM]; /* # of bins in the halo domain */
    int nbinsalloc; /* Total number of allocated bins */
    Bin *bin; /* Array of bins */

    int nstencil; /* # of bins in stencil */
    int *stencil; /* stencil list of bin offsets */

    int iter; /* Iteration of this box */
    
    Domain dombox;  /* Extension of the box in space units */
    Domain domcore; /* Extension of the box minus R_neigh */
    Domain domhalo; /* Extension of the box plus R_neigh */
    Domain dommax;  /* Local atoms must be inside this domain */

    Subdomain sub[NSUB]; /* Array of subdomains */

    int fresh_ghost;  /* 1 it the ghosts have just been recomputed */

    /* These vectors hold the `nlocal` local atoms and, immediately
     * after, the `nghost` ghost atoms. They may not use all the fields
     * (ghost don't have velocity). The allocated size is `nalloc` */
    Vec *r; /* Atom position */
    Vec *v; /* Atom velocity */
    Vec *f; /* Atom force */
    int *atomtype; /* Atom type (mimics original code complexity) */
    Nearby *nearby; /* Nearby atom lists */

    int force_iter; /* Current force iteration */
    int thermo_iter;

    /* Thermal information keep per time step in reduced units */
    double pot_energy; /* Potential energy */
    double kin_energy; /* Kinetic energy */
    double vdwl_energy; /* Van der Waals pairwise energy */
    double virial_pressure; /* Virial pressure */
    double potghost_energy; /* Potential energy of the ghosts only */
    double temperature;
    int ninteractions;

    Hist fhist; /* Force histogram */
    Hist vhist; /* Velocity histogram */
    Hist dhist; /* Nearby atom distance histogram */

    Neigh neigh[NNEIGH]; /* Neighboring boxes info */

    /* Communicators for each type of buffer */
    MPI_Comm comm_r;
    MPI_Comm comm_rt;
    MPI_Comm comm_rvt;

    /* Contiguous pointers to the send and recv PackBuf */
    PackBuf *pb[PB_NTYPES][NNEIGH];

} Box;

typedef struct force {

    /* The R_force^2 for every pair of atom types. Currently all atoms
     * have the same value (R_force^2) but is kept to mimic the
     * complexity of the original code */
    double *R_force_sq;

    /* Parameters for the Lennard-Jones force. Also for each pair of
     * atom types. */
    double *epsilon, *sigma6, *sigma;

} Force;

typedef struct thermo {
    double pot_energy;
} ThermoT;

/* Information for each MPI rank */
typedef struct sim {
    const char *inputfile;
    /* Input data in the order it appears in the input file. */
    //int units;        /* Not used */
    //char *datafile;   /* Not used */
    //char *forcetype;  /* Not used */
    double epsilon, sigma; /* Lennard-Jonnes parameters */
    int npoints[NDIM]; /* Total number or corners of the FCC lattice */
    int timesteps;
    double dt;
    double t_request;   /* Initial temperature */
    double rho;         /* Density */
    int neighbor_period;
    double R_force, skin;
    int thermo_period;
    int nboxes;         /* Boxes per process in Y */
    int nprocsx;
    int nprocsz;

    /* Common constants */
    double mass;
    double dtforce;
    /* 
     * The two radius R_neigh and R_force control the extend of the
     * atomic interaction. The R_force is directly specified in the
     * input file while R_neigh is computed as R_force + skin, where the
     * skin is specified in the input file as well.
     *
     *        . - ^^^ - .        Only the force of the atoms enclosed by
     *     .'     ___     '.     R_force is taken into account. However,
     *    /    `      `     \    all the atoms enclosed by R_neigh are
     *   /   /           \   \   checked, as they could have moved
     *  |   |       R_force   |  inside the R_force region.
     *  |<-------- @ ---->|   |
     *  |  R_neigh        .   |  An example is an atom which moves from
     *   \   \       B   /   /   the point A to the point B, and now
     *    \    -  ___ -     /    is taken into consideration when
     *     `-          A  -'     updating the force of the current atom
     *       '~-. ___ .-~'       (represented by @ at the center).
     *
     */
    double R_neigh;
    double lattice_sep; /* FCC cube side length (equal in all dim) */
    double *R_neigh_sq; /* For each pair of atom types */
    int sort_period;

    /* The potential energy at R_force: this is used to correct the
     * sudden jump in energy when an atom leaves the R_force interaction
     * zone. */
    double e_cut;

    /* Unit conversion constants */
    double mvv2e;
    double dof_boltz;
    double t_scale;
    double p_scale;
    double e_scale;

    int nbinsbox[NDIM]; /* # of bins in the box domain */
    Vec binlen;         /* Length of each bin per dimension */
    Vec invbinlen;      /* Inverse of binlen (to avoid division) */

    int nboxesdim[NDIM];    /* Total number of boxes per dimension */
    int ranknboxesdim[NDIM];    /* # of boxes per dimension of this rank */
    Vec boxlen;
    Vec worldlen;

    /* Process info */
    int rank; /* Current rank index */
    int nranks; /* Total number of ranks */
    int rankdim[NDIM]; /* Cartesian rank coordinates */
    int nranksdim[NDIM]; /* Number of ranks (MPI processes) per dimension */
    MPI_Comm cartesian; /* Cartesian communicator */
    int rankcoord[NDIM]; /* Coordinates of the process */
    
    char *refdir; /* Directory with reference output (NULL disables) */

    double E0_pot; /* Potential energy at start */
    double E0_kin; /* Kinetic energy at start */
    double E0_tot; /* Total energy at start */

    double Epot; /* Current potential energy */
    double Ekin; /* Current kinetic energy */
    double Etot; /* Current total energy */

    int ntypes; /* Number of atom types (species) */
    int ntotatoms;
    int iter;   /* Current iteration from the main task */

    /* GASPI info */
    size_t segsize;
    double *sendseg;
    double *recvseg;
    int nqueues;

    Force force;
    Box *box;
} Sim;

void parse_input(Sim *sim, int argc, char *argv[]);

void box_grow_array(Box *box);
void box_add_atom(Box *box, Vec r, Vec v, int type);

void force_init(Sim *sim);
void force_update(Sim *sim);

void *safe_realloc(void *ptr, size_t size);

void build_nearby_atoms(Sim *sim);

void thermo_update(Sim *sim);
void thermo_init(Sim *sim);

void integrate_init(Sim *sim);
void integrate_position(Sim *sim);
void integrate_velocity(Sim *sim);

#endif /* TYPES_H */
