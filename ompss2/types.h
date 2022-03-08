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

/* Debug checks */

#define MIN_DISTSQ 0.03

/* Extra checks: these allow early catch of problems but introduce large
 * or very large overheads, so they should be disabled when taking into
 * account the performance of the simulation. */

/* Print a histogram of the force magnitudes per box. It should be
 * smooth. */
#define ENABLE_FORCE_HIST
#define FORCE_HIST_NBINS 30
#define FORCE_HIST_MAX 200.0

/* Compute the energy during the simulation. Needed to validate the
 * results. */
#define ENABLE_REALTIME_ENERGY

/* Halts the simulation if the force is too large */
#define ENABLE_MAX_FORCE
#define MAX_FORCE 1e6

/* Halts the simulation if an atom doesn't interact with at least half
 * the neighbors (they are too far away to interact). This may happen
 * with too many time steps without re-neighboring. */
#define ENABLE_MIN_INTERACTIONS

/* Ensure that no new atom is too close to a local atom (slow) */
//#define ENABLE_NEW_ATOM_CHECK

/* Ensure that no ghost atom is too close to a local atom (slow) */
//#define ENABLE_GHOST_ATOM_CHECK


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

typedef struct sim Sim;

typedef struct bin {
    int natoms; /* Number of atoms present in this bin */
    int nalloc;	/* Number of atoms allocated in the bin */
    int *atom; /* Indexes of the atoms */

    /* Thermal information keep per time step in reduced units */
    double *pot_energy; /* Potential energy */
    double *kin_energy; /* Kinetic energy */
    double *vdwl_energy; /* Van der Waals pairwise energy */
    double *virial_temp; /* Virial temperature */
} Bin;

#define NNEIGHSIDE 1
#define NNEIGHDIM (NNEIGHSIDE*2 + 1) 
#define NNEIGH (NNEIGHDIM*NNEIGHDIM*NNEIGHDIM - 1)
#define NSUB 27

typedef struct {
    int natoms;     /* Number of atoms currently in the buffer */
    int nalloc;     /* Number of atoms allocated */
    int atomsize;   /* Number of doubles required per atom */
    double *buf;    /* The contiguous buffer */
    MPI_Request req;
    MPI_Comm comm;
} PackBuf;

typedef struct box Box;
typedef struct neigh Neigh;

/* A neighboring box */
typedef struct neigh {
    int i;              /* Local index of this neighbor in the box */
    int boxcoord[NDIM]; /* Corresponding box coordinate without wrapping */
    int delta[NDIM];    /* Delta vector in boxes */
    int rank;           /* Neighbor process rank */
    int rankcoord[NDIM];/* Neighbor process coordinates without wrapping */

    Box *box;           /* Neighbor box or NULL if outside the rank */
    Neigh *opposite;    /* Opposite neighbor at -delta */

    int wraps;  /* Non-zero if some dimension needs PBC correction */
    Vec addpbc; /* PBC correction per dimension */

    /* Communication packing buffers */

    PackBuf send_rt; /* Send atom positions and types */
    PackBuf recv_rt; /* Receive atom positions and types */

    PackBuf send_rvt; /* Send atom position, velocity and type */
    PackBuf recv_rvt; /* Receive atom positions, velocity and type */
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
    
    Domain dombox;  /* Extension of the box in space units */
    Domain domcore; /* Extension of the box minus R_neigh */
    Domain domhalo; /* Extension of the box plus R_neigh */

    Subdomain sub[NSUB]; /* Array of subdomains */

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

    double *pot_energy;
    double *kin_energy;
    double *virial_temp;

    Neigh neigh[NNEIGH]; /* Neighboring boxes info */
} Box;

typedef struct force {

    /* The R_force^2 for every pair of atom types. Currently all atoms
     * have the same value (R_force^2) but is kept to mimic the
     * complexity of the original code */
    double *R_force_sq;

    /* Parameters for the Lennard-Jones force. Also for each pair of
     * atom types. */
    double *epsilon, *sigma6, *sigma;

    /* Force histogram per box */
    int **forcehist;
    double *forcehistmin;
    double forcehistdelta;

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

    /* Unit conversion constants */
    double mvv2e;
    double dof_boltz;
    double t_scale;
    double p_scale;
    double e_scale;

    int nbinsbox[NDIM];  /* # of bins in the box domain */
    Vec binlen; /* Length of each bin per dimension */

    int nboxesdim[NDIM]; /* Total number of boxes per dimension */
    Vec boxlen;
    Vec worldlen;

    /* Process info */
    int rank; /* Current rank index */
    int nranks; /* Total number of ranks */
    int rankdim[NDIM]; /* Cartesian rank coordinates */
    int nranksdim[NDIM]; /* Number of ranks (MPI processes) per dimension */
    MPI_Comm cartesian; /* Cartesian communicator */
    int rankcoord[NDIM]; /* Coordinates of the process */

    int ntypes; /* Number of atom types (species) */
    int ntotatoms;

    Force force;
    Box *box;
} Sim;

void parse_input(Sim *sim, int argc, char *argv[]);

void box_grow_array(Box *box);
void box_add_atom(Box *box, Vec r, Vec v, int type);

void force_init(Sim *sim);

void comm_setup(Sim *sim);
void comm_atoms_correct_box(Sim *sim);
void comm_borders(Sim *sim);

void *safe_realloc(void *ptr, size_t size);

void packbuf_mpisend(PackBuf *pb, int remoterank, int tag);
void packbuf_mpirecv(PackBuf *pb, int remoterank, int tag);
void packbuf_shmcopy(PackBuf *src, PackBuf *dst);
void packbuf_add_rt(PackBuf *pb, Vec r, int type);
void packbuf_add_rvt(PackBuf *pb, Vec r, Vec v, int type);
void packbuf_unpack_rt(PackBuf *pb, Vec *r, int *types);
void packbuf_unpack_rvt(PackBuf *pb, Vec *r, Vec *v, int *types);
void packbuf_clear(PackBuf *pb);
void packbuf_init(PackBuf *pb, int atomsize);

void build_nearby_atoms(Sim *sim);

#endif /* TYPES_H */
