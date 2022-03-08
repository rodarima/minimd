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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void
force_hist_init(Force *force, int nboxes)
{
#ifdef ENABLE_FORCE_HIST
    force->forcehist = (int **) calloc(nboxes, sizeof(int *));
    force->forcehistmin = (double *) calloc(FORCE_HIST_NBINS, sizeof(double));
    force->forcehistdelta = FORCE_HIST_MAX / (FORCE_HIST_NBINS + 1);

    for (int i = 0; i < nboxes; i++) {
        force->forcehist[i] = (int *) calloc(FORCE_HIST_NBINS, sizeof(int));
    }

    for (int j = 0; j < FORCE_HIST_NBINS; j++) {
        force->forcehistmin[j] = j * force->forcehistdelta;
    }
#endif
}

static void
force_hist_add(Force *force, double *f, int ibox)
{
#ifdef ENABLE_FORCE_HIST
    double forcemag = 0.0;
    for (int d = X; d <= Z; d++) {
        forcemag += f[d] * f[d];
    }
    forcemag = sqrt(forcemag);
    int forcebin = (int) (forcemag / force->forcehistdelta);
    if (forcebin >= FORCE_HIST_NBINS)
        forcebin = FORCE_HIST_NBINS - 1;
    if (forcebin < 0)
        abort();
    force->forcehist[ibox][forcebin]++;
#endif
}

static void
force_hist_clear(Force *force, int ibox)
{
#ifdef ENABLE_FORCE_HIST
    for (int i = 0; i < FORCE_HIST_NBINS; i++) {
        force->forcehist[ibox][i] = 0;
    }
#endif
}

static void
force_hist_print(Force *force, int box_id)
{
#ifdef ENABLE_FORCE_HIST

    /* Only print the histogram of the first box */
    if (box_id != 0)
        return;

    fprintf(stderr, "Force magnitude histogram for box %d:\n", box_id);
    for (int i = 0; i < FORCE_HIST_NBINS; i++) {
        fprintf(stderr, "  [%e .. %e] -> %d atoms\n",
                force->forcehistmin[i],
                force->forcehistmin[i] + force->forcehistdelta,
                force->forcehist[box_id][i]);
    }
#endif
}

static void
check_max_force(double *f)
{
#ifdef ENABLE_MAX_FORCE
    for (int d = X; d <= Z; d++) {
        if (fabs(f[d]) > MAX_FORCE) {
            fprintf(stderr, "force too large: %e %e %e\n", f[X], f[Y], f[Z]);
            abort();
        }
    }
#endif
}

static void
check_min_interactions(int ninteractions, int n)
{
#ifdef ENABLE_INTERACTIONS_CHECK
    if (ninteractions < n / 2) {
        fprintf(stderr, "too few interactions: %d\n", ninteractions);
        abort();
    }
#endif
}

static double
dotprod(double *v, int n)
{
    double sum = 0.0;

    for (int i = 0; i < n; i++) {
        sum += v[i] * v[i];
    }

    return sum;
}

/* Updates the force acting on a given atom at index `i` by taking
 * into account all `n` neighboring atoms in `ineigh` */
static void
update_force_atom(Force *force, Box *box, Bin *bin, int i, int n, int *ineigh, int ntypes)
{
    int type_offset = box->atomtype[i] * ntypes;
    Vec local_f = { 0.0, 0.0, 0.0 };
#ifdef ENABLE_REALTIME_ENERGY
    bin->vdwl_energy = 0.0;
    bin->virial_temp = 0.0;
#endif

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
        if (i == j)
            abort();

        /* Get neighbor atom position */
        Vec rj = {
            box->r[j][X],
            box->r[j][Y],
            box->r[j][Z]
        };

        /* Compute distance vector */
        Vec delta = { ri[X] - rj[X], ri[Y] - rj[Y], ri[Z] - rj[Z] };
        double sqdist = dotprod(delta, 3);
        int type_ij = type_offset + box->atomtype[j];

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

        ninteractions++;

#ifdef ENABLE_REALTIME_ENERGY
        /* Accumulate Van der Waals energy and virial temperature in
         * real time per bin. Energy needs correction to account the
         * R_force approximation. */
        bin->vdwl_energy += (sr6 - 1.0) * sr6eps;
        bin->virial_temp += sqdist * forcemag;
#endif
    }

    double *f = box->f[i];

    for (int d = X; d <= Z; d++) {
        f[d] = local_f[d];
    }

    force_hist_add(force, f, box->i);
    check_max_force(f);
    check_min_interactions(ninteractions, n);
}

/* Update force for all atoms in the given bin index */
static void
update_force_bin(Force *force, Box *box, Bin *bin, int ntypes)
{
    for (int i = 0; i < bin->natoms; i++) {
        /* Compute the actual atom index */
        int iatom = bin->atom[i];

        /* Ignore ghost atoms */
        if (iatom >= box->nlocal)
            continue;

        int *neighs = box->nearby[iatom].atom;
        int numneighs = box->nearby[iatom].natoms;

        update_force_atom(force, box, bin, iatom, numneighs, neighs,
                ntypes);
    }
}

/* Update force for the local atoms in a box */
static void
update_force_box(Force *force, Box *box, int ntypes)
{
    force_hist_clear(force, box->i);

    /* TODO: We may be able to iterate only through the bins in the box
     * domain */
    for (int i = 0; i < box->nbinsalloc; i++) {
        Bin *bin = &box->bin[i];
        update_force_bin(force, box, bin, ntypes);
    }

    /* Increase the iteration for this box */
    box->force_iter++;

    force_hist_print(force, box->i);
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

    force_hist_init(force, sim->nboxes);

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
        update_force_box(&sim->force, box, sim->ntypes);
    }
}
