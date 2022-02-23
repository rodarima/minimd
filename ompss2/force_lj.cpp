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

#include "force_lj.h"
#include "math.h"
#include "stdio.h"
#include <stdlib.h>

#define MAX_FORCE 1e6

#define FORCEHIST_ENABLE
#define FORCEHIST_NBINS 30
#define FORCEHIST_MAX 200.0

ForceLJ::ForceLJ(int ntypes_, int boxes_per_process_)
{
    cutforce = 0.0;
    use_oldcompute = 0;
    reneigh = 1;
    style = FORCELJ;
    ntypes = ntypes_;
    boxes_per_process = boxes_per_process_; // DSM: Multibox change
    bpp = boxes_per_process;

    cutforcesq = new double[ntypes * ntypes];
    epsilon = new double[ntypes * ntypes];
    sigma6 = new double[ntypes * ntypes];
    sigma = new double[ntypes * ntypes];

    forcehist = (int **) calloc(bpp, sizeof(int *));

    forcehistmin = (double *) calloc(FORCEHIST_NBINS, sizeof(double));
    forcehistdelta = FORCEHIST_MAX / (FORCEHIST_NBINS + 1);

    for (int i = 0; i < bpp; i++) {
        forcehist[i] = (int *) calloc(FORCEHIST_NBINS, sizeof(int));
    }

    for (int j = 0; j < FORCEHIST_NBINS; j++) {
        forcehistmin[j] = j * forcehistdelta;
    }

    for (int i = 0; i < ntypes * ntypes; i++) {
        cutforcesq[i] = 0.0;
        epsilon[i] = 1.0;
        sigma6[i] = 1.0;
        sigma[i] = 1.0;
    }

    eng_vdwl = (double *) malloc(boxes_per_process_ * sizeof(double));
    virial = (double *) malloc(boxes_per_process_ * sizeof(double));
    evflag = (int *) malloc(boxes_per_process_ * sizeof(int));
}

ForceLJ::~ForceLJ()
{
    free(eng_vdwl);
    free(virial);
    free(evflag);
}

void ForceLJ::setup()
{
    for (int i = 0; i < ntypes * ntypes; i++)
        cutforcesq[i] = cutforce * cutforce;
}

static double
dotprod(double *v, int n)
{
    double sum = 0.0;

    for (int i=0; i<n; i++) {
        sum += v[i] * v[i];
    }

    return sum;
}

void
forcehist_add(ForceLJ *lj, double *f, int ibox)
{
#ifdef FORCEHIST_ENABLE
    double forcemag = 0.0;
    for (int d = X; d <= Z; d++) {
        forcemag += f[d] * f[d];
    }
    forcemag = sqrt(forcemag);
    int forcebin = (int) (forcemag / lj->forcehistdelta);
    if (forcebin >= FORCEHIST_NBINS)
        forcebin = FORCEHIST_NBINS - 1;
    if (forcebin < 0)
        abort();
    lj->forcehist[ibox][forcebin]++;
#endif
}

void
forcehist_clear(ForceLJ *lj, int ibox)
{
#ifdef FORCEHIST_ENABLE
    for (int i = 0; i < FORCEHIST_NBINS; i++) {
        lj->forcehist[ibox][i] = 0;
    }
#endif
}

void
forcehist_print(ForceLJ *lj, int box_id)
{
#ifdef FORCEHIST_ENABLE
    fprintf(stderr, "Force magnitude histogram for box %d:\n", box_id);
    for (int i = 0; i < FORCEHIST_NBINS; i++) {
        fprintf(stderr, "  [%e .. %e] -> %d atoms\n",
                lj->forcehistmin[i],
                lj->forcehistmin[i] + lj->forcehistdelta,
                lj->forcehist[box_id][i]);
    }
#endif
}

/* Updates the force acting on a given atom at index `i` by taking
 * into account all `n` neighboring atoms in `ineigh` */
void
ForceLJ::update_force_atom(int i, int n, int *ineigh, Atom *atomdata, bool update_energy, int *marker)
{
    int type_offset = atomdata->type[i] * atomdata->ntypes;
    double local_f[3] = { 0.0, 0.0, 0.0 };

    /* Current atom position vector */
    double ri[3] = {
        atomdata->x[i * PAD + X],
        atomdata->x[i * PAD + Y],
        atomdata->x[i * PAD + Z]
    };

    int ninteractions = 0;

    /* This loop is performance critical */
    for (int k = 0; k < n; k++) {
        int j = ineigh[k];

        /* FIXME: the self atom cannot appear in the neighbor list */
        if (i == j)
            abort();

        /* Get neighbor atom position */
        double rj[3] = {
            atomdata->x[j * PAD + X],
            atomdata->x[j * PAD + Y],
            atomdata->x[j * PAD + Z]
        };

        /* Compute distance vector */
        double delta[3] = { ri[X] - rj[X], ri[Y] - rj[Y], ri[Z] - rj[Z] };
        double sqdist = dotprod(delta, 3);
        int type_ij = type_offset + atomdata->type[j];

        /* Ignore far away atoms */
        if (sqdist >= cutforcesq[type_ij])
            continue;

        double sr2 = 1.0 / sqdist;
        double sr6 = sr2 * sr2 * sr2 * sigma6[type_ij];
        double force = 48.0 * sr6 * (sr6 - 0.5) * sr2 * epsilon[type_ij];

        /* Accumulate force for this neighbor */
        local_f[X] += delta[X] * force;
        local_f[Y] += delta[Y] * force;
        local_f[Z] += delta[Z] * force;

        if (atomdata->box_id == 3 && i == 1047 && j == 1352) {
            fprintf(stderr, "XXX pair 1047--1352 force: %e %e %e\n",
                    delta[X] * force,
                    delta[Y] * force,
                    delta[Z] * force);
        }

        for (int d = X; d <= Z; d++) {
            if (fabs(local_f[d]) > MAX_FORCE) {
                fprintf(stderr, "force too large: %e %e %e\n",
                        local_f[X], local_f[Y], local_f[Z]);
                abort();
            }
        }

        ninteractions++;
        marker[i]++;

        if (update_energy) {
//            /* FIXME: Prevents further decomposition within a box */
//            t_eng_vdwl += sr6 * (sr6 - 1.0) * epsilon[type_ij];
//            t_virial += sqdist * force;
        }
    }

    double *f = atomdata->f[i];

    for (int d = X; d <= Z; d++) {
        f[d] = local_f[d];
    }

    forcehist_add(this, f, atomdata->box_id);

    for (int d = X; d <= Z; d++) {
        if (fabs(f[d]) > MAX_FORCE) {
            fprintf(stderr, "force too large: %e %e %e\n", f[X], f[Y], f[Z]);
            abort();
        }
    }

    if (ninteractions < n / 2) {
        fprintf(stderr, "too few interactions: %d\n", ninteractions);
        fprintf(stderr, "atom at xyz=%e %e %e\n", ri[X], ri[Y], ri[Z]);
        abort();
    }

    if (atomdata->box_id == 3 && i == 1047) {
        fprintf(stderr, "XXX atom 1047 force: %e %e %e\n",
                local_f[X], local_f[Y], local_f[Z]);
    }
}

/* Update force for all atoms in the given bin */
void
ForceLJ::update_force_bin(int ibin, Neighbor *nei, Atom *atomdata,
        int *marker)
{
    int n = nei->bincount[ibin];
    int *bins = &nei->bins[ibin * nei->atoms_per_bin];
    //fprintf(stderr, "updating force in bin %d with %d atoms\n", ibin, n);
    if (nei->binchanges)
        abort();

    for (int i = 0; i < n; i++) {
        /* Compute the actual atom index */
        int iatom = bins[i];
        int *neighs = &nei->neighbors[iatom * nei->maxneighs];
        int numneighs = nei->numneigh[iatom];

        /* Ignore ghost atoms */
        if (iatom >= atomdata->nlocal)
            continue;

        update_force_atom(iatom, numneighs, neighs, atomdata, 0, marker);
    }

    if (nei->binchanges)
        abort();
}

/* Update force for all atoms in a box */
void
ForceLJ::update_force_box(Neighbor *nei, Atom *atomdata)
{
    int ntot = atomdata->nlocal + atomdata->nghost;
    int *marker = (int *) calloc(ntot, sizeof(int));

    if (marker == NULL) {
        abort();
    }

    int hist[31] = { 0 };

    forcehist_clear(this, atomdata->box_id);

    for (int i = 0; i < nei->ntotbins; i++) {
        int n = nei->bincount[i];
        if (n > 30)
            n = 30;

        hist[n]++;
        update_force_bin(i, nei, atomdata, marker);
    }

//    fprintf(stderr, "bin occupation for box %d:\n", atomdata->box_id);
//    for (int i = 0; i <= 30; i++) {
//        fprintf(stderr, " %3d: %d\n", i, hist[i]);
//    }
//    fflush(stderr);

    for (int i = 0; i < atomdata->nlocal; i++) {
        if (marker[i] == 0) {
            fprintf(stderr, "force didn't update atom %d\n", i);
            abort();
        }
    }
    for (int i = atomdata->nlocal; i < ntot; i++) {
        if (marker[i] > 0) {
            fprintf(stderr, "force updated ghost %d\n", i);
            abort();
        }
    }

    /* Only print the histogram of the first box */
    if (atomdata->box_id == 0)
        forcehist_print(this, atomdata->box_id);

    free(marker);
}

void
ForceLJ::compute(Atom &atom, Neighbor &nei)
{
    nei.check(atom);
    check_ghost_overlap(&atom);
    update_force_box(&nei, &atom);
}
