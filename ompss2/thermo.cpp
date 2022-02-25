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

#include "thermo.h"
#include "force_lj.h"
#include "integrate.h"
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

Thermo::Thermo() { }
Thermo::~Thermo() { }

void Thermo::setup(double rho_in, Integrate &integrate, Atom &atom, int units)
{
    rho = rho_in;
    ntimes = integrate.ntimes;

    int maxstat;

    if (nstat == 0)
        maxstat = 2;
    else
        maxstat = ntimes / nstat + 2;

    int nboxes = atom.boxes_per_process;

    /* Store the local values per box during the simulation and perform
     * the reductions at the end. Not critical for performance, so we
     * can just use 2 level arrays */

    steparr = (int **) malloc(maxstat * sizeof(int *));
    tmparr  = (double **) malloc(maxstat * sizeof(double *));
    engarr  = (double **) malloc(maxstat * sizeof(double *));
    prsarr  = (double **) malloc(maxstat * sizeof(double *));

    for (int i = 0; i < maxstat; i++) {
        steparr[i] = (int *) malloc(nboxes * sizeof(int));
        tmparr[i]  = (double *) malloc(nboxes * sizeof(double));
        engarr[i]  = (double *) malloc(nboxes * sizeof(double));
        prsarr[i]  = (double *) malloc(nboxes * sizeof(double));
    }

    if (units == LJ) {
        mvv2e = 1.0;
        dof_boltz = (atom.natoms * 3 - 3);
        t_scale = mvv2e / dof_boltz;
        p_scale = 1.0 / 3 / atom.box.xprd / atom.box.yprd / atom.box.zprd;
        e_scale = 0.5;
    } else if (units == METAL) {
        mvv2e = 1.036427e-04;
        dof_boltz = (atom.natoms * 3 - 3) * 8.617343e-05;
        t_scale = mvv2e / dof_boltz;
        p_scale = 1.602176e+06 / 3 / atom.box.xprd / atom.box.yprd / atom.box.zprd;
        e_scale = 524287.985533; // 16.0;
        integrate.dtforce /= mvv2e;
    }
}

void Thermo::compute(int iflag, Atom *atoms[], Force *force, Timer &timer)
{
    double t, eng, p;

    // DSM: nstat is an input file parameter, description: "thermo calculation every this many steps"
    if (iflag > 0 && iflag % nstat)
        return;

    if (iflag == -1 && nstat > 0 && ntimes % nstat == 0)
        return;

    t_act = 0;
    p_act = 0;

    int istep = iflag;

    if (iflag == -1)
        istep = ntimes;

    if (iflag == 0)
        mstat = 0;

    temperature(atoms, mstat);
    energy(atoms, force, mstat);
    pressure(atoms, force, mstat);

    mstat++;

    /* TODO: Print table at the end */
}

/* reduced potential energy */

void Thermo::energy(Atom *atoms[], Force *force, int slot)
{
    int nboxes = atoms[0]->boxes_per_process;

    for (int ib = 0; ib < nboxes; ++ib) {
        #pragma oss task \
            label("Thermo::energy") \
            in(atoms[ib])
            //in(force->eng_vdwl[ib])
        {
            Neighbor &neighbor = *atoms[ib]->neighbor;
            double e_local = 0.0; //force->eng_vdwl[ib];

            if (neighbor.halfneigh) {
                e_local *= 2.0;
            }

            engarr[slot][ib] = e_local * e_scale;
        }
    }
}

/*  reduced temperature */

// DSM Multibox implementation
void Thermo::temperature(Atom *atoms[], int slot)
{
    int nboxes = atoms[0]->boxes_per_process;

    for (int ib = 0; ib < nboxes; ib++) {
        Atom *a = atoms[ib];

        #pragma oss task \
            label("Thermo::temperature") \
            in(a->v)
        {
            double *v = a->v;
            double t_local = 0.0;

            for (int i = 0; i < a->nlocal; i++) {
                double vx = v[i * PAD + 0];
                double vy = v[i * PAD + 1];
                double vz = v[i * PAD + 2];
                
                t_local += (vx * vx + vy * vy + vz * vz) * a->mass;
            }


            tmparr[slot][ib] = t_local * t_scale;
        }
    }
}

/* reduced pressure from virial
   virial = Fi dot Ri summed over own and ghost atoms, since PBC info is
   stored correctly in force array before reverse_communicate is performed */

void Thermo::pressure(Atom *atoms[], Force *force, int slot)
{
    int nboxes = atoms[0]->boxes_per_process;

    /* FIXME: Merge the local collection into the force loop */

    for (int ib = 0; ib < nboxes; ++ib) {
        #pragma oss task \
            label("Thermo::pressure") \
            out(prsarr[slot][ib])
            //in(force->virial[ib])
        {
            Neighbor &neighbor = *atoms[ib]->neighbor;
            //prsarr[slot][ib] = force->virial[ib];
        }
    }
}

/* Used at setup to normalize all velocities to set the temperature of
 * the simulation. Not critical. */

double Thermo::get_global_temperature(Atom *atoms[])
{
    int nboxes = atoms[0]->boxes_per_process;
    double t_local_sum = 0.0;

    for (int ib = 0; ib < nboxes; ib++) {
        Atom *a = atoms[ib];
        double *v = a->v;

//        fprintf(stderr, "velocity=%e %e %e for first atom in box %d\n",
//                    v[0], v[1], v[2], ib);
//        fprintf(stderr, "velocity=%e %e %e for last atom in box %d\n",
//                    v[(a->nlocal-1)*PAD + 0],
//                    v[(a->nlocal-1)*PAD + 1],
//                    v[(a->nlocal-1)*PAD + 2], ib);

        /* The input dependency in(a->v) creates the dependency for
         * &a->v, the address of the pointer a->v, thus it works as a
         * sentinel, which is invariant to relocations or changes in
         * size of a->v. The other tasks that modify the velocity must
         * use out(a->v) as well. */
        #pragma oss task \
            label("Thermo::get_global_temperature() reduction") \
            in(a->v) reduction(+:t_local_sum)
        {
            double t_local = 0.0;

            for (int i = 0; i < a->nlocal; i++) {
                double vx = v[i * PAD + 0];
                double vy = v[i * PAD + 1];
                double vz = v[i * PAD + 2];
                
                t_local += (vx * vx + vy * vy + vz * vz) * a->mass;
            }

            if (isnan(t_local)) {
                fprintf(stderr, "local temp is nan in box %d\n", ib);
                exit(1);
            }

//            fprintf(stderr, "reducing temperature %e for box %d\n",
//                    t_local, ib);

            t_local_sum += t_local;
        }
    }

    /* Wait until the reduction has finished */
    #pragma oss taskwait in(t_local_sum)

    double temp = 0.0;
    MPI_Allreduce(&t_local_sum, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    fprintf(stderr, "reduced temperature is %e\n", temp);

    /* Adjust temperature units */
    temp *= t_scale;

    fprintf(stderr, "corrected temperature is %e\n", temp);

    return temp;
}
