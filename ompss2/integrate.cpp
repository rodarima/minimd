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
#include "dom.h"
#include "hist.h"

#include <math.h>
#include <stdio.h>

/* Ensure the velocity is not too large */
void check_velocity(Sim *sim, Vec v, double dt)
{
    for (int d = X; d <= Z; d++) {
        double dr = v[d] * dt;

        if (dr > sim->boxlen[d]) {
            fprintf(stderr, "atom moving too fast: %e %e %e\n", v[X], v[Y], v[Z]);
            abort();
        }
    }
}

static double
dotprod(Vec v)
{
    double sum = 0.0;

    for (int i = 0; i < NDIM; i++) {
        sum += v[i] * v[i];
    }

    return sum;
}

/* Performs a half-integration updating the velocity and position of the
 * particles of the given box by using the force */
#pragma oss task \
    label("integrate_position_box") \
    in(*(char **)&box->f) \
    inout(*(char **)&box->r) \
    inout(*(char **)&box->v)
static void
integrate_position_box(Sim *sim, Box *box)
{
    for (int i = 0; i < box->nlocal; i++) {
        for (int d = X; d <= Z; d++)
            box->v[i][d] += sim->dtforce * box->f[i][d];

        for (int d = X; d <= Z; d++)
            box->r[i][d] += sim->dt * box->v[i][d];

        if (ENABLE_DOMAIN_CHECK) {
            if (!in_domain(box->r[i], box->dommax)) {
                fprintf(stderr, "box %d: atom %d at %e %e %e moved outside max domain\n",
                        box->i, i, box->r[i][X], box->r[i][Y], box->r[i][Z]);
                abort();
            }
        }

        if (ENABLE_MAX_VELOCITY_CHECK)
            check_velocity(sim, box->v[i], sim->dt);
    }
}

void
integrate_position(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        integrate_position_box(sim, box);
    }
}

/* Finishes the integration by updating the velocity of the local atoms
 * */
#pragma oss task \
    label("integrate_velocity_box") \
    in(*(char **)&box->f) \
    inout(*(char **)&box->v)
static void
integrate_velocity_box(Sim *sim, Box *box)
{
    if (ENABLE_VHIST)
        hist_clear(&box->vhist);

    for (int i = 0; i < box->nlocal; i++) {
        for (int d = X; d <= Z; d++)
            box->v[i][d] += sim->dtforce * box->f[i][d];

        if (ENABLE_VHIST)
            hist_add(&box->vhist, log(1 + dotprod(box->v[i])));
    }

    if (ENABLE_VHIST && box->i == 0)
        hist_print(&box->vhist, sim->iter);
}

void
integrate_velocity(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        integrate_velocity_box(sim, box);
    }
}

void
integrate_init(Sim *sim)
{
    if (ENABLE_VHIST) {
        for (int i = 0; i < sim->nboxes; i++) {
            Box *box = &sim->box[i];
            hist_init(&box->vhist, 200, "vhist.csv", 5.0/200);
        }
    }
}
