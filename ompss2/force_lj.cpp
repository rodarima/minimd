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

ForceLJ::ForceLJ(int ntypes_, int boxes_per_process_)
{
    cutforce = 0.0;
    use_oldcompute = 0;
    reneigh = 1;
    style = FORCELJ;
    ntypes = ntypes_;
    boxes_per_process = boxes_per_process_; // DSM: Multibox change

    cutforcesq = new double[ntypes * ntypes];
    epsilon = new double[ntypes * ntypes];
    sigma6 = new double[ntypes * ntypes];
    sigma = new double[ntypes * ntypes];

    for (int i = 0; i < ntypes * ntypes; i++) {
        cutforcesq[i] = 0.0;
        epsilon[i] = 1.0;
        sigma6[i] = 1.0;
        sigma[i] = 1.0;
    }

    // DSM: Multibox change: one value per box on this process
    eng_vdwl = (double *) malloc(
        boxes_per_process_ * sizeof(double)); // DSM One of the outputs of compute(). Used in energy()
    virial = (double *) malloc(
        boxes_per_process_ * sizeof(double)); // DSM One of the outputs of compute(). Used in pressure()
    evflag = (int *) malloc(boxes_per_process_ * sizeof(int));
}

ForceLJ::~ForceLJ()
{
    // DSM Multibox
    if (eng_vdwl) {
        free(eng_vdwl);
    }
    if (virial) {
        free(virial);
    }
    if (evflag) {
        free(evflag);
    }
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

/* Updates the force acting on a given atom at index `i` by taking
 * into account all `n` neighboring atoms in `ineigh` */
void
ForceLJ::update_force_atom(int i, int n, int *ineigh, Atom *atomdata, bool update_energy)
{
    int type_offset = atomdata->type[i] * atomdata->ntypes;
    double local_f[3] = { 0.0, 0.0, 0.0 };

    /* Current atom position vector */
    double ri[3] = {
        atomdata->x[i * PAD + X],
        atomdata->x[i * PAD + Y],
        atomdata->x[i * PAD + Z]
    };

    /* This loop is performance critical */
    for (int k = 0; k < n; k++) {
        int j = ineigh[k];

        /* Get neighbor atom position */
        double rj[3] = {
            atomdata->x[j * PAD + X],
            atomdata->x[j * PAD + Y],
            atomdata->x[j * PAD + Z]
        };

        /* Compute distance vector */
        double delta[3] = { rj[X] - ri[X], rj[Y] - ri[Y], rj[Z] - ri[Z] };
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

        if (update_energy) {
//            /* FIXME: Prevents further decomposition within a box */
//            t_eng_vdwl += sr6 * (sr6 - 1.0) * epsilon[type_ij];
//            t_virial += sqdist * force;
        }
    }

    double *f = atomdata->f[i];

    f[X] = local_f[X];
    f[Y] = local_f[Y];
    f[Z] = local_f[Z];
}

/* Update force for all atoms in the given bin */
void
ForceLJ::update_force_bin(int ibin, Neighbor *nei, Atom *atomdata)
{
    int natoms = nei->bincount[ibin];
    int *bins = nei->bins;
    for (int i = 0; i < natoms; i++) {
        int *neighs = &nei->neighbors[i * nei->maxneighs];
        int numneighs = nei->numneigh[i];
        update_force_atom(bins[i], numneighs, neighs, atomdata, 0);
    }
}

/* Update force for all atoms in a box */
void
ForceLJ::update_force_box(Neighbor *nei, Atom *atomdata)
{
    for (int i = 0; i < nei->mbins; i++) {
        update_force_bin(i, nei, atomdata);
    }
}

void
ForceLJ::compute(Atom &atom, Neighbor &nei, Comm &comm, int me)
{
    update_force_box(&nei, &atom);
}
