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

#include "force.h"
#include "types.h"
#include "neighbor.h"
#include "hist.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static double
dotprod(Vec v)
{
    double sum = 0.0;

    for (int i = 0; i < NDIM; i++) {
        sum += v[i] * v[i];
    }

    return sum;
}

static void
check_max_force(Vec f)
{
    if (dotprod(f) > MAX_FORCE_SQ) {
        fprintf(stderr, "force too large: %e %e %e\n", f[X], f[Y], f[Z]);
        abort();
    }
}

static void
check_min_interactions(int ninteractions, int n)
{
    if (ninteractions < n / 2) {
        fprintf(stderr, "too few interactions: %d\n", ninteractions);
        abort();
    }
}

/* Updates the force acting on a given atom at index `i` by taking
 * into account all `n` neighboring atoms in `ineigh` */
static void
update_force_atom(Sim *sim, Force *force, Box *box, Bin *bin, int i, int n, int *ineigh, int ntypes)
{
    int type_offset = box->atomtype[i] * ntypes;
    Vec local_f = { 0.0, 0.0, 0.0 };

    /* Current atom position vector */
    Vec ri = {
        box->r[i][X],
        box->r[i][Y],
        box->r[i][Z]
    };

    int ninteractions = 0;

    /* This loop is performance critical */
    for (int k = 0; k < n; k++) {
        int j = ineigh[k];

        /* FIXME: the self atom cannot appear in the neighbor list */
        //if (i == j)
        //    abort();

        /* Get neighbor atom position */
        Vec rj = {
            box->r[j][X],
            box->r[j][Y],
            box->r[j][Z]
        };

        /* Compute distance vector */
        Vec delta = { ri[X] - rj[X], ri[Y] - rj[Y], ri[Z] - rj[Z] };
        double sqdist = dotprod(delta);
        int type_ij = type_offset + box->atomtype[j];

        if (ENABLE_DHIST)
            hist_add(&box->dhist, sqrt(sqdist));

        /* Ignore far away atoms */
        if (sqdist >= force->R_force_sq[type_ij])
            continue;

        double sr2 = 1.0 / sqdist;
        double sr6 = sr2 * sr2 * sr2 * force->sigma6[type_ij];
        double sr6eps = sr6 * force->epsilon[type_ij];
        double forcemag = 48.0 * (sr6 - 0.5) * sr2 * sr6eps;

        /* Accumulate force for this neighbor */
        local_f[X] += delta[X] * forcemag;
        local_f[Y] += delta[Y] * forcemag;
        local_f[Z] += delta[Z] * forcemag;

        if (ENABLE_MIN_INTERACTIONS_CHECK)
            ninteractions++;

        if (ENABLE_REALTIME_ENERGY) {
            /* Accumulate Van der Waals energy and virial temperature in
             * real time per bin. Energy needs correction to account the
             * R_force approximation. */
            double pot = 4.0 * (sr6 - 1.0) * sr6eps;

            if (ENABLE_ECUT_CORRECTION)
                pot -= sim->e_cut;

            bin->vdwl_energy += pot;
            bin->virial_temp += sqdist * forcemag;
        }
    }

    double *f = box->f[i];

    for (int d = X; d <= Z; d++) {
        f[d] = local_f[d];
    }

    if(ENABLE_FHIST)
        hist_add(&box->fhist, log(1 + sqrt(dotprod(f))));

    if (ENABLE_MAX_FORCE_CHECK)
        check_max_force(f);

    if (ENABLE_MIN_INTERACTIONS_CHECK)
        check_min_interactions(ninteractions, n);
}

/* Update force for all atoms in the given bin index */
static void
update_force_bin(Sim *sim, Force *force, Box *box, Bin *bin, int ntypes)
{
    /* Reset energy accumulators per bin */
    if (ENABLE_REALTIME_ENERGY) {
        bin->vdwl_energy = 0.0;
        bin->virial_temp = 0.0;
    }

    for (int i = 0; i < bin->natoms; i++) {
        /* Compute the actual atom index */
        int iatom = bin->atom[i];

        /* Ignore ghost atoms */
        if (iatom >= box->nlocal)
            continue;

        int *neighs = box->nearby[iatom].atom;
        int numneighs = box->nearby[iatom].natoms;

        update_force_atom(sim, force, box, bin, iatom, numneighs, neighs,
                ntypes);
    }
}

static void
dump_atoms(Sim *sim, Box *box)
{
    if (sim->iter == -1) {
        FILE *f = fopen("atompos.csv", "w");
        fprintf(f, "iter,atom,ghost,x,y,z,neigh\n");
        fclose(f);
    }

    FILE *f = fopen("atompos.csv", "a");
    for (int j = 0; j < box->nlocal + box->nghost; j++) {
        Vec *r = &box->r[j];
        int ghost = j >= box->nlocal;
        int nearby = ghost ? 0 : box->nearby[j].natoms;
        fprintf(f, "%d,%d,%d,%e,%e,%e,%d\n",
                sim->iter, j, ghost, (*r)[X], (*r)[Y], (*r)[Z],
                nearby);
    }
    fclose(f);
}

/* Update force for the local atoms in a box */
static void
update_force_box(Sim *sim, Force *force, Box *box, int ntypes)
{
    if(ENABLE_FHIST)
        hist_clear(&box->fhist);

    if(ENABLE_DHIST)
        hist_clear(&box->dhist);

    if (ENABLE_ATOM_TRACKING)
        dump_atoms(sim, box);

    /* TODO: We may be able to iterate only through the bins in the box
     * domain */
    for (int i = 0; i < box->nbinsalloc; i++) {
        Bin *bin = &box->bin[i];
        update_force_bin(sim, force, box, bin, ntypes);
    }

    if(ENABLE_FHIST && box->i == 0)
        hist_print(&box->fhist, sim->iter);

    if(ENABLE_DHIST && box->i == 0)
        hist_print(&box->dhist, sim->iter);

    /* Increase the iteration for this box */
    box->force_iter++;
}

void
force_init(Sim *sim)
{
    Force *force = &sim->force;
    int pairtypes = sim->ntypes * sim->ntypes;

    force->R_force_sq = (double *) calloc(pairtypes, sizeof(double));
    force->epsilon = (double *) calloc(pairtypes, sizeof(double));
    force->sigma6 = (double *) calloc(pairtypes, sizeof(double));

    for (int i = 0; i < pairtypes; i++) {
        /* All values are the same for all types, but we keep them in a
         * vector to mimic the original code complexity */
        force->R_force_sq[i] = sim->R_force * sim->R_force;
        force->epsilon[i] = sim->epsilon;
        force->sigma6[i] = pow(sim->sigma, 6.0);
    }

    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        hist_init(&box->fhist, 100, "fhist.csv", 10.0/100);
        hist_init(&box->dhist, 400, "dhist.csv", 4.0/400.0);
    }

    fprintf(stderr, "force initialized\n");
}

void
force_free(Force *force)
{
    free(force->R_force_sq);
    free(force->epsilon);
    free(force->sigma6);
}

/* Updates the force in all boxes */
void
force_update(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        update_force_box(sim, &sim->force, box, sim->ntypes);
    }
}
