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

#include "types.h"
#include "neigh.h"
#include "hist.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

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
    if (ninteractions < n * 0.5)
        die("too few interactions: %d/%d\n", ninteractions, n);
}

static void
atom_interact(Sim *sim, Box *box, Bin *bin, 
        int type, double sqdist, Vec delta, Vec local_f)
{
    double sr2 = 1.0 / sqdist;
    double sr6 = sr2 * sr2 * sr2 * sim->force.sigma6[type];
    double sr6eps = sr6 * sim->force.epsilon[type];
    double forcemag = 48.0 * (sr6 - 0.5) * sr2 * sr6eps;

    /* Accumulate force for this neighbor */
    for (int d = X; d <= Z; d++)
        local_f[d] += delta[d] * forcemag;

    if (ENABLE_REALTIME_ENERGY) {
        /* Accumulate Van der Waals energy and virial temperature in
         * real time per bin. Energy needs correction to account the
         * R_force approximation. The 4 constant can be factored out, by
         * this is only used for debugging. */
        double pot = 4.0 * (sr6 - 1.0) * sr6eps;

        if (ENABLE_ECUT_CORRECTION)
            pot -= sim->e_cut;

        bin->vdwl_energy += pot;
        bin->virial_pressure += sqdist * forcemag;
    }
}

static void
check_atom_force(Box *box, Vec f, int ninteractions, int nn)
{
    if (ENABLE_FHIST)
        hist_add(&box->fhist, log(1 + sqrt(dotprod(f))));

    if (ENABLE_MAX_FORCE_CHECK)
        check_max_force(f);

    if (ENABLE_MIN_INTERACTIONS_CHECK)
        check_min_interactions(ninteractions, nn);

    if (ENABLE_COUNT_INTERACTIONS) {
        box->ninteractions += ninteractions;

        /* Not sure if we can always guarantee this property, but is
         * useful for debugging */
        if (box->iter == 1 && ninteractions != 54)
            die("expected 54 interactions\n");
    }
}

/* Updates the force acting on a given atom at index `i` by taking
 * into account all nearby atoms within R_force */
static void
update_force_atom(Sim *sim, Box *box, Bin *bin, int i)
{
    Nearby *nearby = &box->nearby[i];
    int type_offset = box->atomtype[i] * sim->ntypes;
    Vec local_f = { 0.0, 0.0, 0.0 };

    /* Current atom position vector */
    Vec ri = { box->r[i][X], box->r[i][Y], box->r[i][Z] };
    int ninteractions = 0;

    /* This loop is _performance critical_ */
    for (int k = 0; k < nearby->natoms; k++) {
        int j = nearby->atom[k];
        int type_ij = type_offset + box->atomtype[j];

        /* Compute distance vector to the nearby atom */
        Vec delta;
        vec_diff(delta, ri, box->r[j]);
        double sqdist = dotprod(delta);

        /* Add to the histogram far away atoms too */
        if (ENABLE_DHIST)
            hist_add(&box->dhist, sqrt(sqdist));

        /* Ignore far away atoms */
        if (sqdist >= sim->force.R_force_sq[type_ij])
            continue;

        atom_interact(sim, box, bin, type_ij, sqdist, delta, local_f);

        ninteractions++;
    }

    /* Set the atom accumulated force */
    double *f = box->f[i];
    for (int d = X; d <= Z; d++)
        f[d] = local_f[d];

    check_atom_force(box, f, ninteractions, nearby->natoms);
}

/* Update force for all atoms in the given bin index */
static void
update_force_bin(Sim *sim, Box *box, Bin *bin)
{
    /* Reset energy accumulators per bin */
    if (ENABLE_REALTIME_ENERGY) {
        bin->vdwl_energy = 0.0;
        bin->virial_pressure = 0.0;
    }

    for (int i = 0; i < bin->natoms; i++) {
        /* Compute the actual atom index */
        int iatom = bin->atom[i];

        /* Ignore ghost atoms */
        if (iatom >= box->nlocal)
            continue;

        update_force_atom(sim, box, bin, iatom);
    }

    /* Accumulate energy from the bin into the box */
    if (ENABLE_REALTIME_ENERGY) {
        box->vdwl_energy += bin->vdwl_energy;
        box->virial_pressure += bin->virial_pressure;
    }
}

/* Updates the force for all atoms in the box, following each in
 * sequential order */
static void
update_force_box_loop(Sim *sim, Box *box)
{
    /* Reset energy accumulators per bin */
    if (ENABLE_REALTIME_ENERGY) {
        for (int i = 0; i < box->nbinsalloc; i++) {
            Bin *bin = &box->bin[i];
            bin->vdwl_energy = 0.0;
            bin->virial_pressure = 0.0;
        }
    }

    /* Ignore ghost atoms, use only nlocal */
    for (int i = 0; i < box->nlocal; i++) {
        Bin *bin = NULL;

        /* The bin is only needed if we want to accumulate energy */
        if (ENABLE_REALTIME_ENERGY) {
            int ibin = get_atom_bin(sim, box, box->r[i]);
            bin = &box->bin[ibin];
        }

        update_force_atom(sim, box, bin, i);
    }

    /* Accumulate bin energy in the box */
    if (ENABLE_REALTIME_ENERGY) {
        for (int i = 0; i < box->nbinsalloc; i++) {
            Bin *bin = &box->bin[i];
            box->vdwl_energy += bin->vdwl_energy;
            box->virial_pressure += bin->virial_pressure;
        }
    }
}

static void
dump_atoms(Sim *sim, Box *box)
{
    #pragma oss taskwait /* for debug */
    if (box->iter == 0) {
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
                box->iter, j, ghost, (*r)[X], (*r)[Y], (*r)[Z],
                nearby);
    }
    fclose(f);

    if (box->iter == 0) {
        FILE *f = fopen("atomneigh.csv", "w");
        fprintf(f, "iter,atom,i,neigh,x,y,z\n");
        fclose(f);
    }

    {
        FILE *f = fopen("atomneigh.csv", "a");
        for (int j = 0; j < box->nlocal; j++) {
            for (int i = 0; i < box->nearby[j].natoms; i++) {
                int neigh = box->nearby[j].atom[i];
                Vec *r = &box->r[neigh];
                fprintf(f, "%d,%d,%d,%d,%e,%e,%e\n",
                        box->iter, j, i, neigh,
                        (*r)[X], (*r)[Y], (*r)[Z]);
            }
        }
        fclose(f);
    }
}

/* Update force for the local atoms in a box */
#pragma oss task \
    label("update_force_box") \
    in(box->iter) \
    in(box->r) \
    in(box->bin) \
    in(box->nearby) \
    inout(box->vdwl_energy, box->virial_pressure) \
    inout(box->f)
static void
update_force_box(Sim *sim, Box *box)
{
    if (ENABLE_FHIST)
        hist_clear(&box->fhist);

    if (ENABLE_DHIST)
        hist_clear(&box->dhist);

    if (ENABLE_ATOM_TRACKING)
        dump_atoms(sim, box);

    /* Clear box energy counters */
    if (ENABLE_REALTIME_ENERGY) {
        box->vdwl_energy = 0.0;
        box->virial_pressure = 0.0;
    }

    box->ninteractions = 0;

    if (ENABLE_FORCE_BY_BINS) {
        for (int i = 0; i < box->nbinsalloc; i++)
            update_force_bin(sim, box, &box->bin[i]);
    } else {
        update_force_box_loop(sim, box);
    }

    if (ENABLE_COUNT_INTERACTIONS)
        fprintf(stderr, "iter %d box %d total interactions %d\n",
                sim->iter, box->i, box->ninteractions);

    if (ENABLE_FHIST && box->i == 0 && sim->rank == 0)
        hist_print(&box->fhist, box->iter);

    if (ENABLE_DHIST && box->i == 0 && sim->rank == 0)
        hist_print(&box->dhist, box->iter);

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
        if (i != 0)
            continue;

        if (ENABLE_FHIST)
            hist_init(&box->fhist, 100, "fhist.csv", 10.0/100);

        if (ENABLE_DHIST)
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
    for (int i = 0; i < sim->nboxes; i++)
        update_force_box(sim, &sim->box[i]);
}
