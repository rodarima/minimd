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


class Neighbor;

class ThreadData {
public:
    ThreadData()
    {
        mpi_me = 0;
        mpi_num_threads = 0;
        omp_me = 0;
        omp_num_threads = 1;
        teams = 1;
    };
    ~ThreadData() {};
    int mpi_me;
    int mpi_num_threads;
    int omp_me;
    int omp_num_threads;
    int teams;
};

class Atom {
public:

    int comm_size, reverse_size, border_size;

    Atom(int ntypes_, int boxes_per_process_);
    ~Atom();
    void addatom(double, double, double, double, double, double);
    void pbc();
    void growarray();

    void copy(int, int);

    void pack_comm(int, int *, double *, int, double, double, double);
    void unpack_comm(int, int, double *);
    void pack_reverse(int, int, double *);
    void unpack_reverse(int, int *, double *);

    int pack_border(int, double *, int *);
    int unpack_border(int, double *);
    int pack_exchange(int, double *);
    int unpack_exchange(int, double *);
    int skip_exchange(double *);

    double *realloc_2d_double_array(double *, int, int, int);
    double *create_2d_double_array(int, int);
    void destroy_2d_double_array(double *);

    int *realloc_1d_int_array(int *, int, int);
    int *create_1d_int_array(int);
    void destroy_1d_int_array(int *);

    void sort(Neighbor &neighbor);
    void sort_ghosts(Neighbor &neighbor);

private:
    int *binpos;
    int *bins;
    double *x_copy;
    double *v_copy;
    int *type_copy;
    int copy_size;
};


#define TIME_TOTAL 0
#define TIME_COMM 1
#define TIME_FORCE 2
#define TIME_NEIGH 3
#define TIME_TEST 4
#define TIME_N 5

class Timer {
public:
    Timer();
    ~Timer();
    void stamp();
    void stamp(int);
    void stamp_extra_start();
    void stamp_extra_stop(int);
    void barrier_start(int);
    void barrier_stop(int);
    double *array;

private:
#ifdef PREC_TIMER
    timespec previous_time, previous_time_extra;
#endif
    double previous_time_d, previous_time_extra_d;
};

class Neighbor {
public:
    int every; // re-neighbor every this often
    int nbinx, nbiny, nbinz; // # of global bins
    double cutneigh; // neighbor cutoff
    double *cutneighsq; // neighbor cutoff squared
    int ncalls; // # of times build has been called
    int max_totalneigh; // largest # of neighbors ever stored

    int *numneigh; // # of neighbors for each atom
    int *neighbors; // array of neighbors of each atom
    int maxneighs; // max number of neighbors per atom
    int halfneigh;

    int ghost_newton;
    int count;
    Neighbor(int ntypes_);
    ~Neighbor();
    int setup(Atom &); // setup bins based on box and cutoff
    void build(Atom &); // create neighbor list
    void check(Atom &);

    Timer *timer;

    ThreadData *threads;

    // Atom is going to call binatoms etc for sorting
    void binatoms(Atom &atom, int count = -1); // bin all atoms

    int *bincount; // ptr to 1st atom in each bin
    int *bins; // ptr to next atom in each bin
    int binchanges;
    int mbins; // total number of bins per box
    int atoms_per_bin;

    /* Number of bins in each dimension */
    int nbins[NDIM];

    /* Number of total bins accounting all dimensions */
    int ntotbins;

    int coord2bin(double, double, double); // mapping atom coord to a bin

private:
    double xprd, yprd, zprd; // box size

    int nmax; // max size of atom arrays in neighbor
    int ntypes; // number of atom types

    int nstencil; // # of bins in stencil
    int *stencil; // stencil list of bin offsets

    int mbinx, mbiny, mbinz;
    int mbinxlo, mbinylo, mbinzlo;
    double binsizex, binsizey, binsizez;
    double bininvx, bininvy, bininvz;

    int resize;

    /* Inner core range in bins. Used to partition the box into two
     * groups of bins: core and shell. So the atoms in the core can
     * begin the force computation without waiting for the exchange of
     * ghosts atoms in the shell */
    Range core_range;

    /* Size of each bin in space units */
    double binlen[NDIM];

    /* Size of each box in space units */
    double boxlen[NDIM];

    double bindist(int, int, int); // distance between binx

    void create_groups();
    void group_bins();
    int find_group(int index[NDIM]);
};

// DSM: Holds pointer to memory buffer of doubles and its maximum size
class AtomBuffer {
public:
    int maxsize; // maximum number of doubles that can be sent by/received into corresponding buffer
    int natoms; // number of atoms within this buffer
    double *buf; // buffer for all comms

    // PBC flags moved here - one set per (send) buffer rather than one
    // per swap Remove separate pbc_flagx, pbc_flagy, pbc_flagz
    // variables and replaced with value to correct by in each
    // dimension. e.g. pbc_x will be one of -box.xprd, 0, box.xprd.
    // Saves having to compute the product every time in packing
    // routines. pbc_any kept TODO: is pbc_any necessary?
    int pbc_any; // whether any PBC on this buffer
    double pbc_x; // PBC correction in x for this swap
    double pbc_y; // same in y
    double pbc_z; // same in z

    // Used for internal buffers only. Required to know where to unpack
    // 1st atom to outside the communicate() function
    int internal_firstrecv;

    void growsend(int); // Increase size of buffer to be able to store at least the given number of items
    void growrecv(int); // As above but different factor for size increase. More suited to recieve buffer pattern

    // Pack this buffer for border communication. Replaces equivalent atom.pack_border()
    int pack_border(int, double, double, double, int, int **, int &);
};

// DSM Multibox: Move majority of Comm attributes into a struct to
// duplicate per box. Better solution would be to rewrite program to use
// multiple Comm object instances but this would require greater
// structural changes.
struct BoxBufs {
    int nswap; // # of swaps to perform

    // DSM These now duplicated per box layer
    int *sendnum[3], *recvnum[3]; // # of atoms to send/recv in each swap
    int *firstrecv[3]; // where to put 1st recv atom in each swap
    int *comm_send_size[3]; // # of values to send in each comm
    int *comm_recv_size[3]; // # of values to recv in each comm
    int *reverse_send_size[3]; // # of values to send in each reverse
    int *reverse_recv_size[3]; // # of values to recv in each reverse

    int *sendproc, *recvproc; // proc to send/recv with at each swap
    int *sendproc_exc, *recvproc_exc; // proc to send/recv with at each swap for safe exchange
    // DSM Multibox change
    int *sendneigh, *recvneigh; // neighbour directions (UP, RIGHT, etc.) to send/recv with at each swap
    MPI_Request request; // outstanding MPI request for this box
    MPI_Comm comm[3]; // communicator for this box, 3 as 1 per communication function (as bufs below)

    // double* buf;
    //  DSM Replaced with 3*24 buffers,
    //        one for each exchange between neighbours
    //    and one for each function that performs communication: communicate(), exchange(), borders()
    //      = 3 functions * 3 layers * 8 neighbours
    //  Indexed by [function][box_layer][neighbour].
    //  3 functions: COMMUNICATE_FUNCTION, EXCHANGE_FUNCTION, BORDERS_FUNCTION
    //  3 box layers: box_id (same layer), box_id-1 (below), box_id+1 (above)
    //  8 neighbours: UP, TOP_RIGHT, RIGHT, etc.

    AtomBuffer bufs_send[3][3][8]; // send buffers for all comm
    AtomBuffer bufs_recv[3][3][8]; // recv buffers for all comm
    int *sendlists[3][8]; // list of atoms for each buffer to send in each swap. Used in borders() and communicate()
    int maxsendlists[3][8]; // size of sendlist

    // DSM buffers for special case: internal memory transfer in y-direction
    // Again 3 buffers, 1 per function that performs communication
    AtomBuffer internal_buf_send_up[3];
    AtomBuffer internal_buf_send_down[3];
    AtomBuffer internal_buf_recv_up[3];
    AtomBuffer internal_buf_recv_down[3];
    int *internal_sendlist_up;
    int internal_maxsendlist_up;
    int *internal_sendlist_down;
    int internal_maxsendlist_down;

    // Used in exchange. Saves recv buffer pointers and sizes from recv
    // tasks to be used in unpack tasks Attributes to avoid going out of
    // scope before tasks are complete 3 box layers, 5 dimensions per
    // swap
    AtomBuffer *buf_recvs1[3][5];
    int nrecvs1[3][5];
    AtomBuffer *buf_recvs2[3][5];
    int nrecvs2[3][5];
    // As above but for borders
    // 3 box layers, 5 dimensions, 2 swaps
    AtomBuffer *buf_recvs_borders[3][5][2];
    int nrecvs_borders[3][5][2];

    // Buffers used exclusively in communicate meighbour tasks
    // implmentation using a single message per neighbour
    AtomBuffer buf_send_single;
    AtomBuffer buf_recv_single;

    // DSM Replaced these with xneg_slab_lo, etc. constants in atom.box struct
    // double* slablo, *slabhi;          // bounds of slabs to send to other procs

    int copy_size;
    int *nsend_thread;
    int *nrecv_thread;
    int *nholes_thread;
    int **exc_sendlist_thread;
    int *send_flag;
    int *maxsend_thread;
    int maxthreads;
    int maxnlocal;
    int nrecv_atoms;
};

class Comm {
public:
    Comm(int); // DSM Multibox: constructor now takes boxes_per_process as an argument
    ~Comm();
    // void free_box_comms(); // DSM Multibox: for freeing box communicators
    int setup(
        double, int, int, Atom *[], int); // DSM Multibox: Array, two process grid int and nonblocking flag arguments
    // Communicate
    void communicate(Atom **); // DSM Multibox: Now takes array of Atom** structs, one per box
    void communicate_blocking(Atom &, int); // DSM Multibox: Blocking implementation. Now also takes int box_id
    void communicate_blocking_isend(Atom &, int);
    void communicate_blocking_alltasks(Atom *, int);
    void communicate_blocking_alltasks_recvfirst(Atom *, int);
    void communicate_nonblocking_neighbourtasks(Atom *, int);
    void communicate_blocking_neighbourtasks(Atom *, int);
    void communicate_blocking_single_message_neighbourtasks(Atom *, int);
    void communicate_nonblocking_alltasks_tampi_iwait(Atom *, int);
    void communicate_nonblocking_neighbourtasks_tampi_iwaitall(Atom *, int);
    void communicate_mpi(Atom *atom, int box_id);
    void communicate_internal(Atom &, int); // communicate() equivalent of exchange_internal()
    void communicate_internal_send(Atom *, AtomBuffer *, AtomBuffer *, int *, int, int);
    // Exchange
    void exchange(Atom **);
    void exchange_nonblocking_neighbourtasks_tampi_iwaitall(Atom *, int);
    void exchange_pack(Atom *, AtomBuffer bufs_send[3][8], AtomBuffer *, AtomBuffer *);
    void exchange_internal(Atom &, int); // DSM Multibox: New function that performs exchange between boxes on same proc
    void exchange_internal_send(AtomBuffer *, AtomBuffer *);
    void exchange_internal_recv(Atom *, AtomBuffer *);
    // Borders
    void borders(Atom **);
    void borders_nonblocking(Atom **);
    void borders_blocking(Atom &, int); // DSM Multibox: Now also takes int box_id
    void borders_blocking_neighbourtasks(Atom *, int);
    void blocking_nonblocking_neighbourtasks_tampi_iwaitall(Atom *, int);
    void borders_internal(Atom &, int); // borders() equivalent of exchange_internal()
    void borders_internal_send(Atom *);
    void borders_pack(Atom *atom);
    // DSM Multibox: changed below buffer reallocation functions to work per box and per buffer
    void growsend(int, int);
    void growrecv(int, int);
    void growlist(int, int, int);
    // Sentinels for tasking
    char *initialIntegrateSentinels;
    char *sortSentinels;
    char *forceComputeSentinels;
    char *finalIntegrateSentinels;
    char *communicateSentinels;
    char *communicateInternalPackSentinels;
    char *communicateInternalUnpackSentinels;
    char **communicatePackSentinels;
    char **communicateSendSentinel;
    char **communicateRecvSentinels;
    char *exchangeSentinels;
    char *exchangePackSentinels;
    char *exchangeInternalSendSentinels;
    char *exchangeInternalRecvSentinels;
    char *exchangePBCSentinels;
    char **exchangeSend1Sentinels;
    char **exchangeRecv1Sentinels;
    char **exchangeSend2Sentinels;
    char **exchangeRecv2Sentinels;
    char *bordersSentinels;
    char *bordersPackSentinels;
    char *bordersUnpackSentinels;
    char *bordersSendSentinels;
    char (*bordersRecvSentinels)[5][2]; // 5 dimensions, 2 swaps per dimension
    char *bordersInternalSentinels;
    char *bordersInternalSendSentinels;
    char *neighbourBuildSentinels;

public:
    int me; // my proc ID
    // DSM: Multibox change.
    int boxes_per_process; // number of boxes on each MPI process
    BoxBufs *boxBufs; // duplicated buffers and counters required per box

    // DSM: Remaining attributes not duplicated in boxBufs, i.e.
    // constants per box (or values relating to broken code)
    // DSM To switch to 26-way communication, procneigh is now 5 dimensional array:
    //   - Sides (original 2 dimensions): x, z
    //   - Corners (2 new dimensions): "top-right" & "bottom-right", "top-left" & "bottom-left"
    //   - Internal communication between boxes: y
    // (Sides+Corners)*3 for box_id+1, box_id, and box_id-1 layers = 4*3
    // = 12 dimensions. NB: Do not need to have
    // separate elements in procneigh for each box layer. Process ranks
    // will be the same on each layer, only target box
    // IDs differ.
    // + 1 for internal communication in y dimension = 13 dimensions total.
    // Exchange in +ve and -ve direction in each dimension => 13*2 =
    // 26-way communication between process neighbours
    int procneigh[5][2]; // my 26 proc neighs (procneigh[1][0:1] will always be my rank as boxes split in y dimension)
    int procgrid[3]; // # of procs in each dim
    int need[5]; // how many procs away needed in each dim
    int boxgrid[3]; // # of boxes in each dim (analogous to procgrid)

    ThreadData *threads; //

    int check_safeexchange; // if sets give warnings if an atom moves further than subdomain size
    int do_safeexchange; // exchange atoms with all subdomains within neighbor cutoff
    Timer *timer;

    // DSM Flag enabling/disabling non-blocking communication mode
    int nonblocking_enabled;
};

class Integrate;

class Thermo {
public:
    int nstat;
    int mstat;
    int ntimes;

    int **steparr;
    double **tmparr;
    double **engarr;
    double **prsarr;

    Thermo();
    ~Thermo();
    void setup(double, Integrate &integrate, Atom &atom, int);
    void temperature(Atom **, int slot);
    double get_global_temperature(Atom *atoms[]);
//    void energy(Atom **, Force *, int slot);
//    void pressure(Atom *atoms[], Force *force, int slot);
//    void compute(int, Atom **, Force *, Timer &);

    double t_act, p_act, e_act;
    double t_scale, e_scale, p_scale, mvv2e, dof_boltz;

    ThreadData *threads;

private:
    double rho;
};

typedef struct sim Sim;

class Integrate {
public:
    double dt;
    double dtforce;
    int ntimes;
    // DSM Multibox: These converted to local variables in the functions where they are used.
    /*int nlocal, nmax;
    double* x, *v, *f, *xold;*/
    double mass;

    int sort_every;

    Integrate();
    ~Integrate();
    void setup();
    // DSM Multibox change: Atom now an array of atoms[] and Neighbor now an attribute of atoms
//    void run(Sim *sim, Atom *atoms[], Force *, Comm &, Thermo &, Timer &);

    ThreadData *threads;
};

typedef struct bin {
    int natoms; /* Number of atoms in this bin */
    int *iatom; /* Indexes of the atoms */

    /* Thermal information keep per time step in reduced units */
    double *pot_energy; /* Potential energy */
    double *kin_energy; /* Kinetic energy */
    double *vdwl_energy; /* Van der Waals pairwise energy */
    double *virial_temp; /* Virial temperature */
} Bin;

#define NNEIGHSIDE 1
#define NNEIGHDIM (NNEIGHSIDE*2 + 1) 
#define NNEIGH (NNEIGHDIM*NNEIGHDIM*NNEIGHDIM - 1)

typedef struct {
    int natoms;     /* Number of atoms currently in the buffer */
    int nalloc;     /* Number of atoms allocated */
    int atomsize;   /* Number of doubles required per atom */
    double *buf;    /* The contiguous buffer */
    MPI_Request req;
    MPI_Comm comm;
} PackBuf;

typedef struct {
    int boxneigh_negative;
    int boxneigh_positive;

} CommT;

typedef struct box Box;
typedef struct neigh Neigh;

/* A neighboring box */
typedef struct neigh {
    int i;              /* Local index of this neighbor in the box */
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

/* All information needed for a box of the simulation */
typedef struct box {
    int i;      /* Box index for this process */
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

    Vec *r; /* Atom positions */
    Vec *v; /* Velocities */
    Vec *f; /* Forces */
    int *atomtype;

    int maxneighs; /* Allocated number of neighbors per atom */
    int *neighbors; /* Neighbor index 2D array of nlocal X maxneighbors */
    int *numneighs; /* Number of neighbors of a given atom */

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

#endif /* TYPES_H */
