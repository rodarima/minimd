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

#define ENABLE_DEBUG 0
#include "log.h"
#include "types.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
    #pragma oss taskwait in(t_local_sum) /* required */

    /* Reduce temperature from all ranks */
    double temp = 0.0;
    MPI_Allreduce(&t_local_sum, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    /* Adjust temperature units */
    temp *= sim->t_scale;

    //fprintf(stderr, "temperature = %e\n", temp);

    return temp;
}

static void
check_final_energy(Sim *sim, double pot, double kin, double tot)
{
    /* Can only check if the two energy corrections are enabled */
    if (ENABLE_ECUT_CORRECTION && ENABLE_NTOTATOMS_CORRECTION) {
        double relerr = fabs((sim->E0_tot - tot) / sim->E0_tot);

        if (relerr >= MAX_ENERGY_REL_ERROR) {
            die("the final total energy diverges too much %e (max %e)\n",
                    relerr, MAX_ENERGY_REL_ERROR);
        }

        err("total energy relative error %e (max %e)\n",
                relerr, MAX_ENERGY_REL_ERROR);
    }
}

#pragma oss task label("thermo_update") \
    in({sim->box[i].iter,               i=0;sim->nboxes}) \
    in({sim->box[i].vdwl_energy,        i=0;sim->nboxes}) \
    in({sim->box[i].virial_pressure,    i=0;sim->nboxes}) \
    in({sim->box[i].temperature,        i=0;sim->nboxes}) \
    in(sim->iter) \
    out(sim->Ekin, sim->Epot, sim->Etot)
static void
thermo_update_internal(Sim *sim, int iter)
{
    double local_vdwl_energy = 0.0;
    double local_virial_pressure = 0.0;
    double local_temperature = 0.0;
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];

        local_vdwl_energy += box->vdwl_energy;
        local_virial_pressure += box->virial_pressure;
        local_temperature += box->temperature;
    }

    double vdwl_energy, virial_pressure, temperature;

    MPI_Reduce(&local_vdwl_energy, &vdwl_energy,
            1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_virial_pressure, &virial_pressure,
            1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_temperature, &temperature,
            1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    /* Correct reduced units */
    vdwl_energy *= sim->e_scale;
    temperature *= sim->t_scale;

    /* Not initialized yet at -1 */
    if (iter == 0) {
        temperature = get_temperature(sim);
    }

    /* Compute energy */
    double pot_energy = vdwl_energy / sim->ntotatoms;
    double kin_energy = temperature * 3.0 / 2.0;
    double tot_energy = pot_energy + kin_energy;

    /* Only the rank 0 prints the report */
    if (sim->rank != 0)
        return;

    sim->Epot = pot_energy;
    sim->Ekin = kin_energy;
    sim->Etot = tot_energy;

    /* Save initial values of energy */
    if (iter == 0) {
        sim->E0_pot = pot_energy;
        sim->E0_kin = kin_energy;
        sim->E0_tot = tot_energy;
    } else if (iter == sim->timesteps) {
        check_final_energy(sim, pot_energy, kin_energy, tot_energy);
    }

    if (ENABLE_REALTIME_ENERGY) {
        FILE *f = fopen("energy.csv", "a");
        fprintf(f, "%d,%e,%e,%e\n", iter, pot_energy, kin_energy, tot_energy);
        fclose(f);
    }
}

void
thermo_update(Sim *sim)
{
    if (ENABLE_REALTIME_ENERGY)
        thermo_update_internal(sim, sim->iter);
}

void
thermo_init(Sim *sim)
{
    /* Zero bin accumulators for energy */
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        /* All the bins, including outside the box */
        for (int j = 0; j < box->nbinsalloc; j++) {
            Bin *bin = &box->bin[j];
            bin->pot_energy = 0.0;
            bin->kin_energy = 0.0;
            bin->vdwl_energy = 0.0;
            bin->virial_pressure = 0.0;
        }
    }

    if (ENABLE_REALTIME_ENERGY) {
        FILE *f = fopen("energy.csv", "w");
        fprintf(f, "iter,Epot,Ekin,Etot\n");
        fclose(f);
    }
}

///* reduced potential energy */
//
//void Thermo::energy(Atom *atoms[], Force *force, int slot)
//{
//    int nboxes = atoms[0]->boxes_per_process;
//
//    for (int ib = 0; ib < nboxes; ++ib) {
//        #pragma oss task \
//            label("Thermo::energy") \
//            in(atoms[ib])
//            //in(force->eng_vdwl[ib])
//        {
//            Neighbor &neighbor = *atoms[ib]->neighbor;
//            double e_local = 0.0; //force->eng_vdwl[ib];
//
//            if (neighbor.halfneigh) {
//                e_local *= 2.0;
//            }
//
//            engarr[slot][ib] = e_local * e_scale;
//        }
//    }
//}
//
///*  reduced temperature */
//
//// DSM Multibox implementation
//void Thermo::temperature(Atom *atoms[], int slot)
//{
//    int nboxes = atoms[0]->boxes_per_process;
//
//    for (int ib = 0; ib < nboxes; ib++) {
//        Atom *a = atoms[ib];
//
//        #pragma oss task \
//            label("Thermo::temperature") \
//            in(a->v)
//        {
//            double *v = a->v;
//            double t_local = 0.0;
//
//            for (int i = 0; i < a->nlocal; i++) {
//                double vx = v[i * PAD + 0];
//                double vy = v[i * PAD + 1];
//                double vz = v[i * PAD + 2];
//                
//                t_local += (vx * vx + vy * vy + vz * vz) * a->mass;
//            }
//
//
//            tmparr[slot][ib] = t_local * t_scale;
//        }
//    }
//}
//
///* reduced pressure from virial
//   virial = Fi dot Ri summed over own and ghost atoms, since PBC info is
//   stored correctly in force array before reverse_communicate is performed */
//
//void Thermo::pressure(Atom *atoms[], Force *force, int slot)
//{
//    int nboxes = atoms[0]->boxes_per_process;
//
//    /* FIXME: Merge the local collection into the force loop */
//
//    for (int ib = 0; ib < nboxes; ++ib) {
//        #pragma oss task \
//            label("Thermo::pressure") \
//            out(prsarr[slot][ib])
//            //in(force->virial[ib])
//        {
//            Neighbor &neighbor = *atoms[ib]->neighbor;
//            //prsarr[slot][ib] = force->virial[ib];
//        }
//    }
//}
//
