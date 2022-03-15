/* ----------------------------------------------------------------------
   miniMD is a simple, parallel molecular dynamics (MD) code.   miniMD is
   an MD microapplication in the Mantevo project at Sandia National
   Laboratories ( http://www.mantevo.org ). The primary
   authors of miniMD are Steve Plimpton (sjplimp@sandia.gov), Paul Crozier
   (pscrozi@sandia.gov) and Christian Trott (crtrott@sandia.gov).

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This library is free software; you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation;
   either version 3 of the License, or (at your option) any later
   version.

   This code is distributed in the hope that it will be useful, but
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

#include "types.h"
#include "neigh.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <fenv.h>
#include <math.h>
#include <float.h>

#ifdef USE_TAMPI
# include <TAMPI.h>
#else
# define MPI_TASK_MULTIPLE MPI_THREAD_MULTIPLE
#endif

static void
setup_ranks(Sim *sim)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &sim->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &sim->nranks);

    /* Rank grid */
    sim->nranksdim[X] = sim->nprocsx;
    sim->nranksdim[Y] = 1;
    sim->nranksdim[Z] = sim->nprocsz;

    /* Create Cartesian MPI process grid */
    int periods[3] = { 1, 1, 1 };
    MPI_Cart_create(MPI_COMM_WORLD, 3, sim->nranksdim, periods, 0, &sim->cartesian);
    MPI_Cart_get(sim->cartesian, 3, sim->nranksdim, periods, sim->rankdim);
}

static void
setup_boxes(Sim *sim)
{
    sim->box = (Box *) calloc(sim->nboxes, sizeof(Box));
    if (sim->box == NULL) {
        perror("calloc failed");
        abort();
    }

    /* Total boxes per dimension */
    sim->nboxesdim[X] = sim->nprocsx;
    sim->nboxesdim[Y] = sim->nboxes;
    sim->nboxesdim[Z] = sim->nprocsz;

    /* Separation between points in the atom lattice */
    sim->lattice_sep = pow((4.0 / sim->rho), (1.0 / 3.0));

    /* Simulation space or world */
    for (int d = X; d <= Z; d++)
        sim->worldlen[d] = sim->npoints[d] * sim->lattice_sep;

    /* Box size is common for all boxes */
    for (int d = X; d <= Z; d++)
        sim->boxlen[d] = sim->worldlen[d] / sim->nboxesdim[d];

    /* Ensure R_neigh doesn't overlap */
    for (int d = X; d <= Z; d++) {
        if (sim->R_neigh * 2 >= sim->boxlen[d]) {
            fprintf(stderr, "error: box too small in dimension %c\n",
                    "XYZ"[d]);
            fprintf(stderr, "R_neigh=%e overlaps with box len=%e\n",
                    sim->R_neigh, sim->boxlen[d]);
            abort();
        }
    }

    /* Setup information per box */
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];

        box->i = i;
        box->idim[X] = sim->rankdim[X];
        box->idim[Y] = i;
        box->idim[Z] = sim->rankdim[Z];

        for (int d = X; d <= Z; d++) {
            /* Box domain */
            box->dombox[d][LO] = sim->boxlen[d] * box->idim[d];
            box->dombox[d][HI] = sim->boxlen[d] * (box->idim[d] + 1);

            /* Add small room for round errors */
            double delta = sim->R_neigh + 1e-6 * sim->boxlen[d];

            /* Box + external halos */
            box->domhalo[d][LO] = box->dombox[d][LO] - delta;
            box->domhalo[d][HI] = box->dombox[d][HI] + delta;

            /* Box - internal halos */
            box->domcore[d][LO] = box->dombox[d][LO] + delta;
            box->domcore[d][HI] = box->dombox[d][HI] - delta;

            /* Box + max distance for local atoms */
            box->dommax[d][LO] = box->dombox[d][LO] - delta * 3.0;
            box->dommax[d][HI] = box->dombox[d][HI] + delta * 3.0;
        }
    }
}

double
get_bindist(int idelta[NDIM], Vec binlen)
{
    Vec del = { 0 };
    double distsq = 0.0;

    for (int d = X; d <= Z; d++) {
        if (idelta[d] > 0)
            del[d] = (idelta[d] - 1) * binlen[d];
        else if (idelta[d] < 0)
            del[d] = (idelta[d] + 1) * binlen[d];

        distsq += del[d] * del[d];
    }

    return distsq;
}

static int
coord2index(int i[NDIM], int n[NDIM])
{
    return i[Z] * n[Y] * n[X] + i[Y] * n[X] + i[X];
}

static void
setup_bin_stencil(Sim *sim, Box *box, int enclosed_nbins[NDIM])
{
    /* Compute and count the actual number of bins in the stencil */
    box->nstencil = 0;

    int i[NDIM];

    /* XXX: Note that we are not taking into account different R_neigh
     * for each pair of atom types. */
    double R_neigh_sq = sim->R_neigh * sim->R_neigh;

    /* Select those bins that fall inside the R_neigh radius */
    for (i[Z] = -enclosed_nbins[Z]; i[Z] <= enclosed_nbins[Z]; i[Z]++) {
        for (i[Y] = -enclosed_nbins[Y]; i[Y] <= enclosed_nbins[Y]; i[Y]++) {
            for (i[X] = -enclosed_nbins[X]; i[X] <= enclosed_nbins[X]; i[X]++) {
                /* Check if the bin is within R_neigh distance */
                if (get_bindist(i, sim->binlen) < R_neigh_sq) {
                    box->stencil[box->nstencil++] = coord2index(i, box->nbinshalo);
                }
            }
        }
    }
    fprintf(stderr, "box = %d nstencil = %d\n", box->i, box->nstencil);
}

static void
setup_bins_box(Sim *sim, Box *box)
{
    int enclosed_nbins[NDIM];
    int total_enclosed = 1;

    box->nbinsalloc = 1;

    for (int d = X; d <= Z; d++) {

        /* Use limits to compute the bins in the halo domain */
        int lo = box->domhalo[d][LO] / sim->binlen[d];
        int hi = box->domhalo[d][HI] / sim->binlen[d];

        /* FIXME: Not sure why this extra bin is needed. It causes the
         * boxes to have different number of halo bins. */
        if (box->domhalo[d][LO] < 0.0)
            lo--;

        /* Extend by 1 to cover stencil */
        lo--;
        hi++;

        /* Compute the number of bins covering the halo domain */
        box->nbinshalo[d] = hi - lo;
        box->nbinsalloc *= box->nbinshalo[d];

        /* The number of bins inside the R_neigh radius */
        enclosed_nbins[d] = sim->R_neigh / sim->binlen[d];

        /* Enlarge if the bin is larger than the radius */
        if (enclosed_nbins[d] * sim->binlen[d] < 1e-6 * sim->R_neigh)
            enclosed_nbins[d]++;

        total_enclosed *= 2 * enclosed_nbins[d] + 1;
    }

    /* Allocate bins and stencil */
    box->bin = (Bin *) malloc(box->nbinsalloc * sizeof(Bin));
    box->stencil = (int *) malloc(total_enclosed * sizeof(int));
    fprintf(stderr, "box = %d total_enclosed = %d\n", box->i,
            total_enclosed);

    for (int i = 0; i < box->nbinsalloc; i++) {
        /* Initialize empty bins */
        Bin *bin = &box->bin[i];
        memset(bin, 0, sizeof(Bin));
    }

    setup_bin_stencil(sim, box, enclosed_nbins);
}

static void
setup_bins(Sim *sim)
{
    /* Fill the number of bins per box */
    double neighscale = 5.0 / 6.0;

    /* Common information for the bins */
    for (int d = X; d <= Z; d++) {
        /* FIXME: npoints the total number of points */
        sim->nbinsbox[d] = neighscale * sim->npoints[d];

        if (sim->nbinsbox[d] == 0) {
            fprintf(stderr, "error: number of bins in dimension %c is 0\n",
                    "XYZ"[d]);
            abort();
        }

        /* Setup bin length */
        sim->binlen[d] = sim->boxlen[d] / (double) sim->nbinsbox[d];
    }

    /* Setup box specific information for the bins */
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        setup_bins_box(sim, box);
    }
}

/* Park/Miller RNG w/out MASKING, so as to be like f90s version */
double
random(int *idum)
{
    int k;
    double ans;

    int IA = 16807;
    int IM = 2147483647;
    double AM = (1.0 / (double) IM);
    int IQ = 127773;
    int IR = 2836;

    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;

    if (*idum < 0)
        *idum += IM;

    ans = AM * (*idum);
    return ans;
}

static void
advance_index(int s[NDIM], int o[NDIM], int subboxdim)
{
    /* Advance index */
    s[X]++;

    /* Wrap subindex */
    for (int d = X; d < Z; d++) {
        if (s[d] == subboxdim) {
            s[d] = 0;
            s[d+1]++;
        }
    }

    /* Propagate wrap to outer index */
    if (s[Z] == subboxdim) {
        s[Z] = 0;
        o[X]++;
    }

    /* Wrap outer index */
    for (int d = X; d < Z; d++) {
        if (o[d] == subboxdim) {
            o[d] = 0;
            o[d+1]++;
        }
    }
}

/* Initialize atoms on FCC lattice */
static void
setup_atoms_box(Sim *sim, Box *box)
{
    box->nlocal = 0;
    box->nghost = 0;
    box->nalloc = 0;

    box->r = NULL;
    box->v = NULL;
    box->f = NULL;
    box->atomtype = NULL;
    box->nearby = NULL;

    /* Determine loop bounds of lattice subsection that overlaps my
     * sub-box insure loop bounds do not exceed nx,ny,nz */

    double dist = 0.5 * sim->lattice_sep;
    if (dist >= sim->R_force) {
        fprintf(stderr, "fatal: atoms are too far away to interact\n");
        exit(1);
    }

    /* Compute the ranges of indexes that may contain an atom */
    Range idom;
    for (int d = X; d <= Z; d++) {
        idom[d][LO] = (box->dombox[d][LO] / dist) - 1;
        idom[d][HI] = (box->dombox[d][HI] / dist) + 1;

        /* Enforce limits */
        idom[d][LO] = MAX(idom[d][LO], 0);
        idom[d][HI] = MIN(idom[d][HI], 2 * sim->npoints[d] - 1);
    }

    /* Each process generates positions and velocities of atoms on FCC
     * sub-lattice that overlaps its box. Only store atoms that fall in
     * my box. Use atom index (generated from lattice coordinates) as
     * unique seed to generate a unique velocity. Exercise RNG between
     * calls to avoid correlations in adjacent atoms */

    double xtmp, ytmp, ztmp, vx, vy, vz;

    int ind[NDIM] = { 0 };
    int s[NDIM] = { 0 };
    int o[NDIM] = { 0 };
    int subboxdim = 8;

    /* TODO: This initialization can probably be simplified even further
     * to three nested loops iterating only in the finer lattice */

    for (; o[Z] * subboxdim <= idom[Z][HI]; advance_index(s, o, subboxdim)) {
        /* Compute current index vector */
        for (int d = X; d <= Z; d++)
            ind[d] = o[d] * subboxdim + s[d];

        /* Place an atom only in the positions of the face centered
         * cubic lattice */
        if ((ind[X] + ind[Y] + ind[Z]) % 2 != 0)
            continue;

        int skip = 0;

        /* Check boundaries */
        for (int d = X; d <= Z; d++) {
            if (ind[d] < idom[d][LO] || ind[d] > idom[d][HI] + 2) {
                skip = 1;
                break;
            }
        }

        if (skip)
            continue;

        /* Compute atom position in space units */
        Vec r;
        for (int d = X; d <= Z; d++)
            r[d] = dist * ind[d];

        /* Ensure the atom falls inside the box */
        for (int d = X; d <= Z; d++) {
            if (r[d] < box->dombox[d][LO] || r[d] >= box->dombox[d][HI]) {
                skip = 1;
                break;
            }
        }

        if (skip) {
            continue;
        }


        /* Compute deterministic seed for pseudorandom velocity */
        int seed = ind[Z] * (2 * sim->npoints[Y]) * (2 * sim->npoints[X])
            + ind[Y] * (2 * sim->npoints[Y])
            + ind[X] + 1;

        Vec v;
        for (int d = X; d <= Z; d++) {
            for (int m = 0; m < 5; m++)
                random(&seed);

            v[d] = random(&seed);
        }

        int type = rand() % sim->ntypes;

        /* Place atom here */
        box_add_atom(box, r, v, type);
    }
}

static void
check_natoms(Sim *sim)
{
    if (!ENABLE_ATOM_COUNT_CHECK)
        return;

    #pragma oss taskwait
    /* Ensure the total number of atoms is correct */
    int global_natoms = 0;
    int rank_natoms = 0;
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        rank_natoms += box->nlocal;
    }

    MPI_Reduce((void *) &rank_natoms, (void *) &global_natoms, 1,
            MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (sim->rank == 0) {
        if (global_natoms != sim->ntotatoms) {
            fprintf(stderr, "error: total atoms %d mismatch, expected %d\n",
                    global_natoms, sim->ntotatoms);
            abort();
        }
    }
}

static void
setup_atoms(Sim *sim)
{
    /* Setup information per box */
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        fprintf(stderr, "setting atoms for box %d\n", i);
        setup_atoms_box(sim, box);
    }

    check_natoms(sim);
}

static double
get_temperature(Sim *sim)
{
    double t_local_sum = 0.0;

    for (int ib = 0; ib < sim->nboxes; ib++) {
        Box *box = &sim->box[ib];

        /* The input dependency in(box->v) creates the dependency for
         * &box->v, the address of the pointer box->v, thus it works as
         * a sentinel, which is invariant to relocations or changes in
         * size of box->v. The other tasks that modify the velocity must
         * use out(box->v) as well. */
        #pragma oss task \
            label("get_global_temperature() reduction") \
            in(box->v) reduction(+:t_local_sum)
        {
            double t_local = 0.0;

            for (int i = 0; i < box->nlocal; i++) {
                double vx = box->v[i][X];
                double vy = box->v[i][Y];
                double vz = box->v[i][Z];
                
                t_local += (vx * vx + vy * vy + vz * vz) * sim->mass;
            }

            if (isnan(t_local)) {
                fprintf(stderr, "local temp is nan in box %d\n", ib);
                abort();
            }

            t_local_sum += t_local;
        }
    }

    /* Wait until the reduction has finished */
    #pragma oss taskwait in(t_local_sum)

    /* Reduce temperature from all ranks */
    double temp = 0.0;
    MPI_Allreduce(&t_local_sum, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    /* Adjust temperature units */
    temp *= sim->t_scale;

    fprintf(stderr, "temperature = %e\n", temp);

    return temp;
}

/* Adjust initial velocities to give desired temperature */

static void
setup_temperature(Sim *sim)
{
    /* Zero center-of-mass motion */
    Vec vlocal = { 0 };
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < box->nlocal; j++) {
            for (int d = X; d <= Z; d++) {
                vlocal[d] += box->v[j][d];
            }
        }
    }

    if (isnan(vlocal[X] + vlocal[Y] + vlocal[Z])) {
        fprintf(stderr, "local sum of velocities is nan\n");
        abort();
    }

    Vec vtot, vmean;
    for (int d = X; d <= Z; d++) {
        MPI_Allreduce(&vlocal[d], &vtot[d], 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
        vmean[d] = vtot[d] / sim->ntotatoms;
    }

    fprintf(stderr, "vmean = (%e %e %e)\n", vmean[X], vmean[Y], vmean[Z]);

    /* Zero the mean velocity by shifting each atom vmean */
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < box->nlocal; j++) {
            for (int d = X; d <= Z; d++) {
                box->v[j][d] -= vmean[d];
            }
        }
    }

    /* Adjust the temperature by scaling the atoms velocity */
    double t = get_temperature(sim);
    double factor = sqrt(sim->t_request / t);
//    /* FIXME: remove me */
//    factor = 10.0;

    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < box->nlocal; j++) {
            for (int d = X; d <= Z; d++) {
                box->v[j][d] *= factor;
            }
        }
    }

    /* Ensure the temperature is now correct */
    double t_corrected = get_temperature(sim);
    double relerr = fabs(t_corrected - sim->t_request) / fabs(sim->t_request);

    /* This holds when relerr is nan too */
    if (! (relerr < 10e2 * DBL_EPSILON)) {
        fprintf(stderr, "temperature relative error %e (t_corrected=%e vs t_requested=%e)\n",
                relerr, t_corrected, sim->t_request);
        abort();
    }
}

static void
setup_params(Sim *sim)
{
    /* Enable floating point exceptions */
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    /* Set a fixed known seed */
    srand(5413);

    /* Set the current iteration to -1, before sim_run() */
    sim->iter = -1;

    /* Derived constants */
    sim->mass = 1.0; /* Default mass. TODO check optimized value */
    sim->dtforce = 0.5 * sim->dt / sim->mass;
    sim->R_neigh = sim->R_force + sim->skin;
    sim->ntypes = 4;

    /* Given the total number of points in the grid N the total number
     * of atoms can be computed as 4*N. This is a face centered cubic
     * lattice (FCC).
     *
     *     *---------*
     *    /    *    /|
     *   * --------* |
     *   |         |*|
     *   |    *    | *
     *   |         |/
     *   *---------*
     */
    /* FIXME: remove me */
    sim->ntotatoms = 4 * sim->npoints[X] * sim->npoints[Y] * sim->npoints[Z];
    //sim->ntotatoms = 2;

    /* Unit conversion constants */
    sim->mvv2e = 1.0;
    sim->dof_boltz = (sim->ntotatoms * 3 - 3);
    sim->t_scale = sim->mvv2e / sim->dof_boltz;
    sim->p_scale = 1.0 / 3.0 / sim->boxlen[X] / sim->boxlen[Y] / sim->boxlen[Z];
    sim->e_scale = 0.5;

    /* TODO: We should use an array to mimic the complexity of multiple
     * atom types */
    sim->e_cut = 4 * (pow(1.0 / sim->R_force, 12.0)
            - pow(1.0 / sim->R_force, 6.0));

    sim->sort_period = sim->neighbor_period;

    int pairtypes = sim->ntypes * sim->ntypes;

    sim->R_neigh_sq = (double *) calloc(pairtypes, sizeof(double));

    for (int i = 0; i < pairtypes; i++) {
        sim->R_neigh_sq[i] = sim->R_neigh * sim->R_neigh;
    }
}

static void
setup_neighbors_box(Sim *sim, Box *box)
{
    int nn = NNEIGHSIDE;
    int delta[NDIM];

    /* Most of the neighbor information is identical in every box, but
     * is more clear to place it inside each box, along with the buffers
     * for communication. So a Neigh is mapped to the neighbor box. */

    /* Setup neighbor deltas by iterating through all the combinations;
     * it is easier to compute than the other way around. */
    for (delta[Z] = -nn; delta[Z] <= nn; delta[Z]++) {
        for (delta[Y] = -nn; delta[Y] <= nn; delta[Y]++) {
            for (delta[X] = -nn; delta[X] <= nn; delta[X]++) {
                int in = delta2neigh(delta);
                /* Skip the current box */
                if (in < 0)
                    continue;

                Neigh *neigh = &box->neigh[in];
                neigh->i = in;
                for (int d = X; d <= Z; d++) {
                    neigh->delta[d] = delta[d];
                }

                /* Get neighbor rank and box coordinates */
                for (int d = X; d <= Z; d++) {
                    neigh->boxcoord[d] = box->idim[d] + delta[d];
                    if (neigh->boxcoord[d] < 0) {
                        neigh->rankcoord[d] = -1;
                        neigh->boxcoordw[d] = sim->nboxesdim[d] - 1;
                    } else if (neigh->boxcoord[d] >= sim->nboxesdim[d]) {
                        neigh->rankcoord[d] = sim->rankdim[d];
                        neigh->boxcoordw[d] = 0;
                    } else {
                        neigh->rankcoord[d] = neigh->boxcoord[d] / sim->nboxesdim[d];
                        neigh->boxcoordw[d] = neigh->boxcoord[d];
                    }
                }

                /* Use the rank coordinates to find the rank */
                MPI_Cart_rank(sim->cartesian, neigh->rankcoord, &neigh->rank);
            }
        }
    }

    for (int in = 0; in < NNEIGH; in++) {
        Neigh *neigh = &box->neigh[in];
        fprintf(stderr, "neigh %2d (%2d) rank %2d delta (%2d %2d %2d)\n",
                in, neigh->i, neigh->rank, 
                neigh->delta[X], neigh->delta[Y], neigh->delta[Z]);
    }

    /* Set wrapping */
    for (int in = 0; in < NNEIGH; in++) {
        Neigh *neigh = &box->neigh[in];

        /* Initially assume the neighbor doesn't need wrapping */
        neigh->wraps = 0;

        for (int d = X; d <= Z; d++) {
            if (neigh->boxcoord[d] < 0) {
                neigh->addpbc[d] = +sim->worldlen[d];
                neigh->wraps = 1;
            } else if (neigh->boxcoord[d] >= sim->nboxesdim[d]) {
                neigh->addpbc[d] = -sim->worldlen[d];
                neigh->wraps = 1;
            } else {
                neigh->addpbc[d] = 0.0;
            }
        }
    }

    /* Add pointer to opposite neighbor and boxes if they are in the
     * same rank */
    for (int in = 0; in < NNEIGH; in++) {
        Neigh *neigh = &box->neigh[in];

        neigh->opposite = &box->neigh[opposite_neigh(in)];

        if (neigh->rank == sim->rank) {
            /* Only the Y dimension is relevant */
            int ib = box->i + neigh->delta[Y];
            if (ib < 0)
                ib += sim->nboxes;
            else if (ib >= sim->nboxes)
                ib -= sim->nboxes;

            neigh->box = &sim->box[ib];
        } else {
            neigh->box = NULL;
        }
    }
}

static void
setup_neighbors(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++)
        setup_neighbors_box(sim, &sim->box[i]);
}

static int
subdomain_has_neigh(Sim *sim, Subdomain *sub, Neigh *neigh)
{
    for (int d = X; d <= Z; d++) {
        /* Select neighbors in the enclosed volume */
        int lo = MIN(sub->delta[d], 0);
        int hi = MAX(sub->delta[d], 0);

        /* Ignore if out of range */
        if (neigh->delta[d] < lo || neigh->delta[d] > hi)
            return 0;

//        /* Ignore neighbors with single rank dimensions */
//        if (neigh->delta[d] != 0 && sim->nboxesdim[d] == 1)
//            return 0;
    }

    return 1;
}

static void
setup_subdomains_box(Sim *sim, Box *box)
{
    int delta[NDIM];

//    fprintf(stderr, "setup subdomains for box %d\n", box->i);

    for (delta[Z] = -1; delta[Z] <= 1; delta[Z]++) {
        for (delta[Y] = -1; delta[Y] <= 1; delta[Y]++) {
            for (delta[X] = -1; delta[X] <= 1; delta[X]++) {
                int i = delta2subdom(delta);

                Subdomain *sub = &box->sub[i];
                sub->i = i;
                for (int d = X; d <= Z; d++)
                    sub->delta[d] = delta[d];

                sub->nneigh = 0;

                for (int j = 0; j < NNEIGH; j++) {
                    Neigh *neigh = &box->neigh[j];

                    if (!subdomain_has_neigh(sim, sub, neigh))
                        continue;

//                    /* Ensure the neighbor is not already in the list;
//                     * this may happen if by wrapping we end up in the
//                     * same box. */
//                    int skip = 0;
//                    int *coord = neigh->boxcoordw;
//                    for (int k = 0; k < sub->nneigh; k++) {
//                        /* Assume coordinates are equal */
//                        int same = 1;
//                        for (int d = X; d <= Z; d++) {
//                            if (coord[d] != sub->neigh[k]->boxcoordw[d])
//                                same = 0;
//                        }
//
//                        /* If we found a neighbor with the same box
//                         * coordinates, skip it */
//                        if (same) {
//                            skip = 1;
//                            break;
//                        }
//                    }
//
//                    if (skip)
//                        continue;

                    sub->neigh[sub->nneigh++] = neigh;

//                    fprintf(stderr, "box %d sub (%2d %2d %2d): adding neigh (%2d %2d %2d)\n",
//                            box->i,
//                            sub->delta[X], sub->delta[Y], sub->delta[Z],
//                            neigh->delta[X], neigh->delta[Y], neigh->delta[Z]);
                }
            }
        }
    }
}

static void
setup_subdomains(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++)
        setup_subdomains_box(sim, &sim->box[i]);
}

static void
setup_packbuf(Sim *sim)
{
    /* Init all pack buffers */
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < NNEIGH; j++) {
            Neigh *neigh = &box->neigh[j];

            /* Setup the number of doubles needed per buffer */
            packbuf_init(&neigh->send_r,    0, NDIM);
            packbuf_init(&neigh->recv_r,    0, NDIM);
            packbuf_init(&neigh->send_rt,   1, NDIM + 1);
            packbuf_init(&neigh->recv_rt,   1, NDIM + 1);
            packbuf_init(&neigh->send_rvt,  0, NDIM + NDIM + 1);
            packbuf_init(&neigh->recv_rvt,  0, NDIM + NDIM + 1);

            packbuf_debug_switch(&neigh->send_r, PB_GARBAGE, PB_READY);
            packbuf_debug_switch(&neigh->recv_r, PB_GARBAGE, PB_READY);
            packbuf_debug_switch(&neigh->send_rt, PB_GARBAGE, PB_READY);
            packbuf_debug_switch(&neigh->recv_rt, PB_GARBAGE, PB_READY);
            packbuf_debug_switch(&neigh->send_rvt, PB_GARBAGE, PB_READY);
            packbuf_debug_switch(&neigh->recv_rvt, PB_GARBAGE, PB_READY);
        }
    }
}

static void
print_params(Sim *sim)
{
    if (sim->rank != 0)
        return;

    /* Print the input file and details, so the output becomes a valid
     * input file */

    FILE *f = fopen(sim->inputfile, "r");

    if (f == NULL) {
        perror("fopen failed");
        abort();
    }

    while (!feof(f)) {
        char line[4096];
        size_t nread = fread(line, 1, sizeof(line), f);

        printf("nread = %ld\n", nread);

        if (nread == 0) {
            perror("fread failed");
            abort();
        }

        if (fwrite(line, 1, nread, stdout) != nread) {
            perror("fwrite failed");
            abort();
        }
    }
    fclose(f);

    printf("#\n");
    printf("# %s output ...\n", VARIANT_STRING);
    printf("#\n");
    printf("# Run Settings: \n");
    printf("#   MPI processes: %d\n", sim->nranks);
    printf("#   Inputfile: %s\n", sim->inputfile);
    printf("#   Boxes per process: %d\n", sim->nboxes); // DSM Multibox change
    printf("#   Process Grid: %d x %d x %d\n",
            sim->nranksdim[X], sim->nranksdim[Y], sim->nranksdim[Z]);
    printf("# Physics Settings: \n");
    printf("#   ForceStyle: LJ\n");
    printf("#   Force Parameters: epsilon=%2.2lf sigma=%2.2lf\n", sim->epsilon, sim->sigma);
    printf("#   Units: LJ\n");
    printf("#   Atoms: %d\n", sim->ntotatoms);
    printf("#   Atom types: %d\n", sim->ntypes);
    printf("#   System size: (x=%2.2lf y=%2.2lf z=%2.2lf)\n",
            sim->worldlen[X], sim->worldlen[Y], sim->worldlen[Z]);
    printf("#   Unit cells: (x=%i y=%i z=%i)\n",
            sim->npoints[X], sim->npoints[Y], sim->npoints[Z]);
    printf("#   Density: %lf\n", sim->rho);
    printf("#   Force cutoff radius: %lf\n", sim->R_force);
    printf("#   Timestep size: %lf\n", sim->dt);
    printf("# Technical Settings: \n");
    printf("#   Neighbor cutoff radius: %lf\n", sim->R_neigh);
    printf("#   Neighbor bins in box domain: (%i %i %i)\n",
            sim->nbinsbox[X], sim->nbinsbox[Y], sim->nbinsbox[Z]);
    printf("#   Re-neighbor period: %i\n", sim->neighbor_period);
    printf("#   Sorting period: %i\n", sim->sort_period);
    printf("#   Thermo period: %i\n", sim->thermo_period);

    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        fprintf(stderr, "box %d has %d local atoms\n", box->i, box->nlocal);
    }
}

void
sim_init(Sim *sim, int argc, char *argv[])
{
    /* Read input file and set input values */
    parse_input(sim, argc, argv);

    /* Compute derived constants from the input */
    setup_params(sim);

    /* MPI process related parameters */
    setup_ranks(sim);

    /* Create boxes */
    setup_boxes(sim);

    /* Prepare bins sizes and parameters */
    setup_bins(sim);

    /* Create atoms and place then in the lattice with random velocity */
    setup_atoms(sim);

    /* Scale the atom velocity to match the requested initial temperature */
    setup_temperature(sim);

    /* Setup force parameters */
    force_init(sim);

    /* Setup integrator */
    integrate_init(sim);

    setup_neighbors(sim);

    setup_subdomains(sim);

    /* Init pack buffers */
    setup_packbuf(sim);

    thermo_init(sim);

    print_params(sim);

    /* Move atoms to their correct box.
     * FIXME: this should be unneeded, as the atoms must be already
     * initialized in their correct box. */
    comm_tidy(sim);
    #pragma oss taskwait

    /* Copy the ghost atoms into the neighbor processes */
    comm_borders(sim);
    #pragma oss taskwait

    build_nearby_atoms(sim);
    #pragma oss taskwait

    force_update(sim);
    #pragma oss taskwait

    thermo_update(sim);
    #pragma oss taskwait
}

void
sim_run(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        box->iter = 0;
    }

    /* Main simulation loop */
    for (sim->iter = 0; sim->iter < sim->timesteps; sim->iter++) {
        fprintf(stderr, "===== RUNNING ITERATION %d =====\n", sim->iter);

        int recompute_neigh = ((sim->iter + 1) % sim->neighbor_period == 0);
        int print_thermo_stats = ((sim->iter + 1) % sim->thermo_period == 0);

        /* Update atoms positions and half velocities */
        integrate_position(sim);
        check_natoms(sim);

        if (!recompute_neigh) {
            comm_ghost_position(sim);
        } else {
            comm_tidy(sim);
            check_natoms(sim);
            //sort_atoms(sim);
            comm_borders(sim);
            check_natoms(sim);
            build_nearby_atoms(sim);
        }

//    #pragma oss taskwait /* Fixes the problem */
        check_natoms(sim);
        force_update(sim);
        integrate_velocity(sim);

        if (print_thermo_stats)
            thermo_update(sim);

        for (int i = 0; i < sim->nboxes; i++) {
            Box *box = &sim->box[i];
            /* Wait for the velocity or thermo to finish */
            #pragma oss task in(box->v) inout(box->iter)
            box->iter++;
        }
    }

    #pragma oss taskwait
}

void
sim_finalize()
{
//    // DSM Multibox: Sum nlocal over all boxes in this process and contribute that to the Allreduce.
//    // DSM TODO: What is the purpose of this Allreduce? natoms appears to only be used in the PERF_SUMMARY
//    int natoms;
//    int all_boxes_nlocal = 0;
//    for (int box_index = 0; box_index < in.nboxes; ++box_index) {
//        all_boxes_nlocal += atoms[box_index]->nlocal;
//    }
//    MPI_Allreduce(&all_boxes_nlocal, &natoms, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//
//    force_update(&sim, atoms);
//
//    // DSM TODO: Changed to get code to compile - this path will be broken for multibox
//    if (atoms[0]->neighbor->halfneigh && atoms[0]->neighbor->ghost_newton) {
//        fprintf(stderr, "halfneigh and ghost_newton not supported\n");
//        MPI_Finalize();
//        exit(1);
//    }
//
//    thermo.compute(-1, atoms, NULL, timer);
//
//    if (me == 0) {
//        double time_other
//            = timer.array[TIME_TOTAL] - timer.array[TIME_FORCE] - timer.array[TIME_NEIGH] - timer.array[TIME_COMM];
//        printf("\n\n");
//        printf("# Performance Summary:\n");
//        printf("# MPI_proc OMP_threads nsteps natoms t_total t_force t_neigh t_comm t_other performance perf/thread "
//               "grep_string t_extra\n");
//        printf("%i %i %i %i %lf %lf %lf %lf %lf %lf %lf PERF_SUMMARY %lf\n\n\n", nprocs, num_threads, integrate.ntimes,
//            natoms, timer.array[TIME_TOTAL], timer.array[TIME_FORCE], timer.array[TIME_NEIGH], timer.array[TIME_COMM],
//            time_other, 1.0 * natoms * integrate.ntimes / timer.array[TIME_TOTAL],
//            1.0 * natoms * integrate.ntimes / timer.array[TIME_TOTAL] / nprocs / num_threads, timer.array[TIME_TEST]);
//    }
}

int main(int argc, char *argv[])
{
    /* Don't zero so we can find error with asan */
    Sim *sim = (Sim *) malloc(sizeof(Sim));
    int provided;

    MPI_Init_thread(&argc, &argv, MPI_TASK_MULTIPLE, &provided);

    if (provided != MPI_TASK_MULTIPLE) {
        fprintf(stderr, "error: MPI_TASK_MULTIPLE not supported\n");
        abort();
    }

    /* Initialize the simulation structures using the input
     * configuration */
    sim_init(sim, argc, argv);
    sim_run(sim);
    //sim_finalize(sim);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
