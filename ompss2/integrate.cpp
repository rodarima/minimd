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
//#define PRINTDEBUG(a) a
#define PRINTDEBUG(a)
#include "integrate.h"
#include "math.h"
#include <stdio.h>

Integrate::Integrate() { sort_every = 20; }
Integrate::~Integrate() { }

void Integrate::setup() { dtforce = 0.5 * dt; }

/* Ensure the position is within a reasonable limit */
void check_position(Vec r, Atom *atom)
{
    double factor = 2.0;
    Box *box = &atom->box;

    for (int d=X; d<=Z; d++) {
        double lo = box->dom[d][LO] - factor * box->len[d];
        double hi = box->dom[d][HI] + factor * box->len[d];

        if (r[d] < lo || r[d] > hi) {
            fprintf(stderr, "box %d: atom too far: %e %e %e\n",
                    atom->box_id, r[X], r[Y], r[Z]);
            abort();
        }

//        if (r[d] < box->dom[d][LO] || r[d] > box->dom[d][HI]) {
//            fprintf(stderr, "warning: box %d, atom out of box: %e %e %e\n",
//                    atom->box_id, r[X], r[Y], r[Z]);
//        }
    }
}

/* Ensure the velocity is not too large */
void check_velocity(Vec v, double dt, Atom *atom)
{
    Box *box = &atom->box;

    for (int d=X; d<=Z; d++) {
        double dr = v[d] * dt;

        if (dr > box->len[d]) {
            fprintf(stderr, "atom moving too fast: %e %e %e\n", v[X], v[Y], v[Z]);
            abort();
        }
    }
}

/* Performs a half-integration updating the velocity and position of the
 * particles of the given box by using the force */
void initial_integrate(Atom *atoms[], double dt, double dtforce)
{
    int nboxes = atoms[0]->boxes_per_process;
    for (int ib = 0; ib < nboxes; ib++) {
        Atom* a = atoms[ib];

        /* Set the dependencies over the pointers f, v and x (not the
         * region towards they point to) so they serve as sentinels */
        #pragma oss task label("initial_integrate") \
            firstprivate(a) in(a->f) inout(a->v) inout(a->x)
        {
            double *x = a->x;
            double *v = a->v;
            double (*f)[PAD] = a->f;
            size_t n = a->nlocal;
            size_t pad = PAD;

            fprintf(stderr, "initial_integrate for box %d\n", ib);

            for (int i = 0; i < n; i++) {

                v[i * PAD + 0] += dtforce * f[i][X];
                v[i * PAD + 1] += dtforce * f[i][Y];
                v[i * PAD + 2] += dtforce * f[i][Z];

                x[i * PAD + 0] += dt * v[i * PAD + 0];
                x[i * PAD + 1] += dt * v[i * PAD + 1];
                x[i * PAD + 2] += dt * v[i * PAD + 2];

                if (i == 1047 && ib == 3) {
                    fprintf(stderr, "XXX initial integrate v = %e %e %e\n",
                            v[i * PAD + 0],
                            v[i * PAD + 1],
                            v[i * PAD + 2]);
                    fprintf(stderr, "XXX initial integrate f = %e %e %e\n",
                            f[i][X],
                            f[i][Y],
                            f[i][Z]);
                }

                check_position(&x[i * PAD + 0], a);
                check_velocity(&v[i * PAD + 0], dt, a);
            }
        }
    }
}

/* Finishes the integration by updating the velocity of all particles */
void final_integrate(Atom *atoms[], double dtforce)
{
    int nboxes = atoms[0]->boxes_per_process;
    for (int ib = 0; ib < nboxes; ib++) {

        Atom* a = atoms[ib];

        /* Set the dependencies over the pointers f and v (not the
         * region towards they point to) so they serve as sentinels */
        #pragma oss task label("final_integrate for half velocity") \
            firstprivate(a) in(a->f) inout(a->v)
        {
            double *x = a->x;
            double *v = a->v;
            double (*f)[PAD] = a->f;
            size_t n = a->nlocal;
            size_t pad = PAD;
            fprintf(stderr, "final_integrate for box %d\n", ib);

            for (int i = 0; i < n; i++) {
                v[i * PAD + 0] += dtforce * f[i][X];
                v[i * PAD + 1] += dtforce * f[i][Y];
                v[i * PAD + 2] += dtforce * f[i][Z];
            }
        }
    }
}

void sort_atoms(Atom *atoms[], Comm *comm)
{
    int nboxes = atoms[0]->boxes_per_process;

    for (int i = 0; i < nboxes; i++) {
        Atom* a = atoms[i];

        #pragma oss task \
            label("atom->sort") \
            in(comm->exchangePBCSentinels[i]) \
            out(comm->sortSentinels[i]) \
            firstprivate(a)
        {
            a->sort(*a->neighbor);
        }
    }
}

void neigh_build(Atom *atoms[], Comm *comm)
{
    int nboxes = atoms[0]->boxes_per_process;

    for (int i = 0; i < nboxes; i++) {
        Atom* a = atoms[i];

        // Depends on all borders unpack tasks completing,
        // i.e. must have full knowledge of ghost atoms
        // before rebuilding neighbour list.
        // Does not modify x, no pack depedencies.
        // Dependencies on communicate tasks to prevent this
        // task running before all non-rebuild iterations
        // are complete
        #pragma oss task \
            label("neighbor->build")                                        \
            in(comm->bordersUnpackSentinels[i]) \
            in(comm->bordersInternalSentinels[i]) \
            in(comm->communicateSentinels[i]) \
            in(comm->communicateInternalUnpackSentinels[i]) \
            out(comm->forceComputeSentinels[i]) \
            firstprivate(a)
        {
            a->neighbor->build(*a);
        }
    }
}

void force_compute(Atom *atoms[], Comm *comm, Force *force, int print_thermo_stats)
{
    int nboxes = atoms[0]->boxes_per_process;

    for (int i = 0; i < nboxes; i++) {
        Atom* a = atoms[i];
        // No need for borders dependencies, borders tasks run only in
        // reneighbouring branch
        #pragma oss task \
            label("force->compute") \
            in(comm->communicateSentinels[i]) \
            in(comm->communicateInternalUnpackSentinels[i]) \
            out(comm->forceComputeSentinels[i]) \
            out(a->f) \
            out(force->eng_vdwl[i]) \
            out(force->virial[i]) \
            firstprivate(i, a)
        {
            // DSM: thermo.nstat is a constant, an input file parameter fixed at initial setup
            force->evflag[i] = print_thermo_stats;
            // Controls whether eng_vdwl and virial are set this
            // compute call or not.
            // The last 2 arguments (comm & comm.me) are not used in
            // force_lj implementation. Replace with nulls
            force->compute(*a, *a->neighbor);
        }
    }
}

void Integrate::run(Atom *atoms[], Force *force, Comm &comm, Thermo &thermo, Timer &timer)
{
    comm.timer = &timer;
    timer.array[TIME_TEST] = 0.0;

    int check_safeexchange = comm.check_safeexchange;
    const int every
        = (*atoms[0]->neighbor)
              .every; // DSM Multibox: "every" is constant across all neighbors and we can assume at least 1 box
    const int boxes_per_process = atoms[0]->boxes_per_process;

    char *initialIntegrateSentinels = comm.initialIntegrateSentinels;
    char *sortSentinels = comm.sortSentinels;
    char *exchangePBCSentinels = comm.exchangePBCSentinels;
    char *communicateSentinels = comm.communicateSentinels;
    char **communicatePackSentinels = comm.communicatePackSentinels;
    char *communicateInternalPackSentinels = comm.communicateInternalPackSentinels;
    char *communicateInternalUnpackSentinels = comm.communicateInternalUnpackSentinels;
    char *bordersInternalSentinels = comm.bordersInternalSentinels;
    char *forceComputeSentinels = comm.forceComputeSentinels;
    char *neighbourBuildSentinels = comm.neighbourBuildSentinels;
    char *bordersUnpackSentinels = comm.bordersUnpackSentinels;
    char *exchangePackSentinels = comm.exchangePackSentinels;
    char *bordersPackSentinels = comm.bordersPackSentinels;

    // DSM Multibox change: assuming all atoms have the same mass value here and we have at least 1
    mass = atoms[0]->mass;
    dtforce = dtforce / mass;

    // Replaced with an array (one per box)
    int next_sort[atoms[0]->boxes_per_process];
    for (int i = 0; i < atoms[0]->boxes_per_process; ++i) {
        next_sort[i] = sort_every > 0 ? sort_every : ntimes + 1;
    }

    #pragma oss taskwait

    for (int iter = 0; iter < ntimes; iter++) {
        int recompute_neigh = ((iter + 1) % every == 0);
        int print_thermo_stats = ((iter + 1) % thermo.nstat == 0);

        /* Update atoms positions and half velocities */
        initial_integrate(atoms, dt, dtforce);

        if (!recompute_neigh) {
            #pragma oss taskwait
            comm.communicate(atoms);
            #pragma oss taskwait
        } else {
            /* expensive */
            #pragma oss taskwait
            comm.exchange(atoms);
            #pragma oss taskwait
            sort_atoms(atoms, &comm);
            #pragma oss taskwait
            comm.borders(atoms);
            #pragma oss taskwait
            neigh_build(atoms, &comm);
            #pragma oss taskwait
        }

        #pragma oss taskwait
        force_compute(atoms, &comm, force, print_thermo_stats);
        #pragma oss taskwait
        final_integrate(atoms, dtforce);

        if (print_thermo_stats)
            thermo.compute(iter + 1, atoms, force, timer);
    }

    #pragma oss taskwait
}
