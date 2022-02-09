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

#include "atom.h"
#include "integrate.h"
#include "mpi.h"
#include "neighbor.h"
#include "thermo.h"
#include "types.h"
#include <float.h>
#include <cmath>
#include <cstdio>

#include <cstdio>
#include <cstring>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

double random(int *);

#define NSECTIONS 3
#define MAXLINE 255
char line[MAXLINE];
char keyword[MAXLINE];
FILE *fp;

/* create simulation box */

void create_box(Atom &atom, int nx, int ny, int nz, double rho)
{
    double lattice = pow((4.0 / rho), (1.0 / 3.0));
    atom.box.xprd = nx * lattice;
    atom.box.yprd = ny * lattice;
    atom.box.zprd = nz * lattice;
}

/* initialize atoms on fcc lattice in parallel fashion */

int create_atoms(Atom &atom, int nx, int ny, int nz, double rho)
{
    /* total # of atoms */
    // DSM: Why 4*? Is this related to crystal latice structure?
    // RAM: There are 4 species of atoms (types)
    atom.natoms = 4 * nx * ny * nz;
    atom.nlocal = 0;

    /* determine loop bounds of lattice subsection that overlaps my sub-box
       insure loop bounds do not exceed nx,ny,nz */

    double alat = pow((4.0 / rho), (1.0 / 3.0));
    int ilo = static_cast<int>(atom.box.xlo / (0.5 * alat) - 1);
    int ihi = static_cast<int>(atom.box.xhi / (0.5 * alat) + 1);
    int jlo = static_cast<int>(atom.box.ylo / (0.5 * alat) - 1);
    int jhi = static_cast<int>(atom.box.yhi / (0.5 * alat) + 1);
    int klo = static_cast<int>(atom.box.zlo / (0.5 * alat) - 1);
    int khi = static_cast<int>(atom.box.zhi / (0.5 * alat) + 1);

    ilo = MAX(ilo, 0);
    ihi = MIN(ihi, 2 * nx - 1);
    jlo = MAX(jlo, 0);
    jhi = MIN(jhi, 2 * ny - 1);
    klo = MAX(klo, 0);
    khi = MIN(khi, 2 * nz - 1);

    /* each proc generates positions and velocities of atoms on fcc sublattice
         that overlaps its box
       only store atoms that fall in my box
       use atom # (generated from lattice coords) as unique seed to generate a
         unique velocity
       exercise RNG between calls to avoid correlations in adjacent atoms */

    double xtmp, ytmp, ztmp, vx, vy, vz;
    int i, j, k, m, n;
    int sx = 0;
    int sy = 0;
    int sz = 0;
    int ox = 0;
    int oy = 0;
    int oz = 0;
    int subboxdim = 8;

    int iflag = 0;

    while (oz * subboxdim <= khi) {
        k = oz * subboxdim + sz;
        j = oy * subboxdim + sy;
        i = ox * subboxdim + sx;

        if (iflag)
            continue;

        if (((i + j + k) % 2 == 0) && (i >= ilo) && (i <= ihi) && (j >= jlo) && (j <= jhi) && (k >= klo)
            && (k <= khi)) {

            // DSM: Atom position calculations
            xtmp = 0.5 * alat * i;
            ytmp = 0.5 * alat * j;
            ztmp = 0.5 * alat * k;

            // DSM: "only store atoms that fall in my box" = no need to coordinate atom ownership via MPI.
            if (xtmp >= atom.box.xlo && xtmp < atom.box.xhi && ytmp >= atom.box.ylo && ytmp < atom.box.yhi
                && ztmp >= atom.box.zlo && ztmp < atom.box.zhi) {
                n = k * (2 * ny) * (2 * nx) + j * (2 * nx) + i + 1; // DSM: This is the random seed/atom #

                // DSM: Velocity (v) calculations
                for (m = 0; m < 5; m++)
                    random(&n);

                vx = random(&n);

                for (m = 0; m < 5; m++)
                    random(&n);

                vy = random(&n);

                for (m = 0; m < 5; m++)
                    random(&n);

                vz = random(&n);

                atom.addatom(xtmp, ytmp, ztmp, vx, vy, vz);
            }
        }

        sx++;

        if (sx == subboxdim) {
            sx = 0;
            sy++;
        }

        if (sy == subboxdim) {
            sy = 0;
            sz++;
        }

        if (sz == subboxdim) {
            sz = 0;
            ox++;
        }

        if (ox * subboxdim > ihi) {
            ox = 0;
            oy++;
        }

        if (oy * subboxdim > jhi) {
            oy = 0;
            oz++;
        }
    }

    /* check for overflows on any proc */
    // DSM TODO: does this work as intended? int iflag = 0; is before while loop above and only ever read, not written
    // All ranks will always contribute iflag=0 to the allreduce - redundant call?
    // Is this a check for memory corruption with the assumption that iflag will probably be overwritten if we have an
    // overflow? Multibox change: removing these checks as they're not strictly necessary and would not be
    // straightforward to fix. Could have each process do a local sum/max over boxes and contribute that to the
    // Allreduces.
    /*
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    int iflagall;
    MPI_Allreduce(&iflag, &iflagall, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    if(iflagall) {
      if(me == 0) printf("No memory for atoms\n");

      return 1;
    }

    // check that correct # of atoms were created
    // DSM: atom.nlocal incremented by each atom.addatom() call
    int natoms;
    MPI_Allreduce(&atom.nlocal, &natoms, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(natoms != atom.natoms) {
      if(me == 0) printf("Created incorrect # of atoms\n");

      return 1;
    }*/

    return 0;
}

/* adjust initial velocities to give desired temperature */

void create_velocity(double t_request, Atom *atoms[], Thermo &thermo)
{
    int i;

    /* zero center-of-mass motion */

    double vxtot = 0.0;
    double vytot = 0.0;
    double vztot = 0.0;

    /* FIXME: Use tasks? */

    // DSM Multibox change: Perform local sum over all boxes
    for (int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
        for (i = 0; i < atoms[box_index]->nlocal; i++) {
            vxtot += atoms[box_index]->v[i * PAD + 0];
            vytot += atoms[box_index]->v[i * PAD + 1];
            vztot += atoms[box_index]->v[i * PAD + 2];
        }
    }

    if (isnan(vxtot + vytot + vztot)) {
        fprintf(stderr, "local sum of velocities is nan\n");
        exit(1);
    }

    double tmp;
    MPI_Allreduce(&vxtot, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    vxtot = tmp / atoms[0]->natoms; // natoms is constant over all boxes
    MPI_Allreduce(&vytot, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    vytot = tmp / atoms[0]->natoms;
    MPI_Allreduce(&vztot, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    vztot = tmp / atoms[0]->natoms;

    fprintf(stderr, "v/natoms = (%e %e %e)\n",
            vxtot, vytot, vztot);

    // DSM Multibox change
    for (int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
        for (i = 0; i < atoms[box_index]->nlocal; i++) {
            atoms[box_index]->v[i * PAD + 0] -= vxtot;
            atoms[box_index]->v[i * PAD + 1] -= vytot;
            atoms[box_index]->v[i * PAD + 2] -= vztot;
        }
    }

    /* rescale velocities, including old ones */
    thermo.t_act = 0;
    double t = thermo.get_global_temperature(atoms);
    double factor = sqrt(t_request / t);

    fprintf(stderr, "correcting temperature from %e to %e\n",
            t, t_request);

    // DSM Multibox change
    for (int box_index = 0; box_index < atoms[0]->boxes_per_process; ++box_index) {
        for (i = 0; i < atoms[box_index]->nlocal; i++) {
            atoms[box_index]->v[i * PAD + 0] *= factor;
            atoms[box_index]->v[i * PAD + 1] *= factor;
            atoms[box_index]->v[i * PAD + 2] *= factor;
        }
    }

    /* Ensure the temperature is now correct */
    double t2 = thermo.get_global_temperature(atoms);
    fprintf(stderr, "corrected temperature from %e to %e (requested %e)\n",
            t, t2, t_request);

    double relerr = fabs(t2 - t_request) / fabs(t_request);
    fprintf(stderr, "initial temperature relative error %e\n", relerr);

    /* This holds when relerr is nan too */
    if (! (relerr < 10 * DBL_EPSILON)) {
        fprintf(stderr, "temperature relative error %e (t=%e vs treq=%e)\n",
                relerr, t2, t_request);
        exit(1);
    }
}

/* Park/Miller RNG w/out MASKING, so as to be like f90s version */

#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double random(int *idum)
{
    int k;
    double ans;

    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;

    if (*idum < 0)
        *idum += IM;

    ans = AM * (*idum);
    return ans;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK
