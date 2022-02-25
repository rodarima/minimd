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

typedef double Vec[NDIM];
typedef double Domain[NDIM][NLIM];
typedef int    Range[NDIM][NLIM];

#define MIN_DISTSQ 0.03

typedef struct bin {
    /* Thermal information keep per time step */
    double *pot_energy;
    double *kin_energy;
    double *virial_temp;
} Bin;


/* All information needed for a box of the simulation */
typedef struct box {

    /* Information about the bins */
    int nbins;
    Bin *bin;

    /* Bins per dimension */
    int nbinsdim[NDIM];

    /* Box dimensions in space units */
    Vec len;

    /* Extension of the box with the center being 0 */
    Domain dom;

    /* -- Housekeeping -- */

    int force_iter; /* Current force iteration */

    /* -- Kinematics -- */

    //Vec *r; /* Atom positions */
    //Vec *v; /* Velocities */
    //Vec *f; /* Forces */

    /* FIXME: Replace all these names with Ranges */

    /* Legacy information */
    double xprd, yprd, zprd;
    double xlo, xhi;
    double ylo, yhi;
    double zlo, zhi;

    // DSM 26-way communication change: Replace comm's slablo/slabhi
    // arrays with these constants (set in comm.setup()) The
    // comm.borders() function uses these to determine whether atoms are
    // in the slab region before updating positions. Corner calculations
    // use constants from multiple dimensions, e.g. communication with
    // "top-left" neighbour process is all atoms in region boxed by
    // (xneg_slab_lo to xneg_slab_hi) and (zneg_slab_lo to
    // zneg_slab_hi).
    double xneg_slab_lo, xneg_slab_hi, xpos_slab_lo, xpos_slab_hi;
    double yneg_slab_lo, yneg_slab_hi, ypos_slab_lo, ypos_slab_hi;
    double zneg_slab_lo, zneg_slab_hi, zpos_slab_lo, zpos_slab_hi;

} Box;


typedef struct input {
    /* Input data as it appears in the input file */
    int nx, ny, nz;
    double t_request;
    double rho;
    int units;
    ForceStyle forcetype;
    double epsilon, sigma;
    char *datafile;
    int ntimes;
    double dt;
    int neigh_every;
    double R_force;
    double skin_len;
    int thermo_nstat;
    int boxes_per_process;
    int nprocsx;
    int nprocsz;
    int nonblocking_enabled;

    /* Computed / constant */
    double R_neigh;
    int ntypes; /* Number of atom types (species) */
} Input;

typedef struct force {

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

/* Information for each MPI rank */
typedef struct sim {
    Input input;
    Force force;
} Sim;

extern Sim sim;

#endif /* TYPES_H */
