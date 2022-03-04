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
#include "types.h"
#include "mpi.h"
#include "neighbor.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <math.h>

#define DELTA 20000

void check_new_atom(Atom *a, int iatom, Vec ri) {
#ifdef ENABLE_NEW_ATOM_CHECK

    /* Compare with all local atoms */

    for (int j = 0; j < a->nlocal; j++) {
        /* Cannot be */
        if (iatom == j)
            abort();

        double *rj = &a->x[j * PAD];
        double distsq = 0.0;
        for (int d = X; d <= Z; d++) {
            distsq += (ri[d] - rj[d]) * (ri[d] - rj[d]);
        }

        if (distsq < MIN_DISTSQ) {
            fprintf(stderr, "box %d: new atom %d(%d) too close to local atom %d\n",
                    a->box_id, iatom, iatom - a->nlocal, j);
            abort();
        }
    }
#endif
}

/* Ensures that ghost atoms are far from local atoms in O(n^2) */
void check_ghost_overlap(Atom *a)
{
#ifdef ENABLE_GHOST_ATOM_CHECK
    /* Compare with all local atoms */
    for (int i = 0; i < a->nlocal; i++) {
        double *ri = &a->x[i * PAD];
        for (int j = a->nlocal; j < a->nlocal + a->nghost; j++) {
            /* Cannot be */
            if (i == j)
                abort();

            double *rj = &a->x[j * PAD];
            double distsq = 0.0;

            for (int d = X; d <= Z; d++) {
                distsq += (ri[d] - rj[d]) * (ri[d] - rj[d]);
            }

            if (distsq < MIN_DISTSQ) {
                fprintf(stderr, "box %d: local atom %d too close to ghost atom %d\n",
                        a->box_id, i, j);
                abort();
            }
        }
    }
#endif
}

void Atom::growarray()
{
    int nold = nmax;
    nmax += DELTA;
    x = (double *) realloc_2d_double_array(x, nmax, PAD, PAD * nold);
    v = (double *) realloc_2d_double_array(v, nmax, PAD, PAD * nold);
    f = (double (*)[PAD]) realloc_2d_double_array((double *) f, nmax, PAD, PAD * nold);
    type = realloc_1d_int_array(type, nmax, nold);
    xold = (double *) realloc_2d_double_array(xold, nmax, PAD, PAD * nold);

    if (x == NULL || v == NULL || f == NULL || xold == NULL) {
        printf("ERROR: No memory for atoms\n");
    }
}

// DSM: Called in loop at initial setup from setup.cpp's create_atoms()
// DSM: x_in, y_in, z_in are atom coordinates. vx_in, vy_in, vz_in are atom velocities
void Atom::addatom(double x_in, double y_in, double z_in, double vx_in, double vy_in, double vz_in)
{
    // DSM: If we're out of memory for new atoms, allocate arrays that are "DELTA" elements larger, copy existing
    // contents over and free the originals. On first call, nlocal and nmax are 0 - initial allocation is done without
    // any copies.
    if (nlocal == nmax)
        growarray();

    // DSM: x-array is atom position, v-array is atom velocity
    x[nlocal * PAD + 0] = x_in;
    x[nlocal * PAD + 1] = y_in;
    x[nlocal * PAD + 2] = z_in;
    v[nlocal * PAD + 0] = vx_in;
    v[nlocal * PAD + 1] = vy_in;
    v[nlocal * PAD + 2] = vz_in;
    type[nlocal] = rand() % ntypes; // DSM: Counts for each atom type randomly decided here, not in input file.

    nlocal++;
}

/* enforce PBC
   order of 2 tests is important to insure lo-bound <= coord < hi-bound
   even with round-off errors where (coord +/- epsilon) +/- period = bound */

void Atom::pbc()
{
    //#pragma omp for
    for (int i = 0; i < nlocal; i++) {
        if (x[i * PAD + 0] < 0.0)
            x[i * PAD + 0] += box.xprd;

        if (x[i * PAD + 0] >= box.xprd)
            x[i * PAD + 0] -= box.xprd;

        if (x[i * PAD + 1] < 0.0)
            x[i * PAD + 1] += box.yprd;

        if (x[i * PAD + 1] >= box.yprd)
            x[i * PAD + 1] -= box.yprd;

        if (x[i * PAD + 2] < 0.0)
            x[i * PAD + 2] += box.zprd;

        if (x[i * PAD + 2] >= box.zprd)
            x[i * PAD + 2] -= box.zprd;
    }
}

void Atom::copy(int i, int j)
{
    x[j * PAD + 0] = x[i * PAD + 0];
    x[j * PAD + 1] = x[i * PAD + 1];
    x[j * PAD + 2] = x[i * PAD + 2];
    v[j * PAD + 0] = v[i * PAD + 0];
    v[j * PAD + 1] = v[i * PAD + 1];
    v[j * PAD + 2] = v[i * PAD + 2];
    type[j] = type[i];
}

// DSM: Packs box's position array x into buf for use in MPI_Send. PBC corrections applied here.
// DSM: Parameters:
//   - n          number of atoms/size of list
//   - list       array of indices/atom IDs into the x/atom position array. These are the atoms to be sent.
//   - buf        buffer atoms are packed into
//   - pbc_flags  periodic boundary conditions array
void Atom::pack_comm(int n, int *list, double *buf, int pbc_any, double pbc_x, double pbc_y, double pbc_z)
{
    int i, j;

    // DSM: Why does this branch exist? pbc_flags[1-3] can be -1, 1 or 0 depending on whether the box length is to be
    // subtracted, added or not included in the calculation (i.e. atom is within the box). Top and bottom branches give
    // identical answers when pbc_flags[1-3] == 0, i.e. when pbc_flags[0] == 0.
    if (pbc_any == 0) {
        for (i = 0; i < n; i++) {
            j = list[i];
            double *r = &x[j * PAD];
//            fprintf(stderr, "box %d: packing atom %d at %e %e %e to send\n",
//                    this->box_id, j, r[X], r[Y], r[Z]);
            buf[3 * i] = x[j * PAD + 0];
            buf[3 * i + 1] = x[j * PAD + 1];
            buf[3 * i + 2] = x[j * PAD + 2];
        }
    } else {
        for (i = 0; i < n; i++) {
            j = list[i];
            double *r = &x[j * PAD];
//            fprintf(stderr, "box %d: packing atom %d at %d pbc old pos %e %e %e\n",
//                    this->box_id, j, i, r[X], r[Y], r[Z]);
            buf[3 * i] = x[j * PAD + 0] + pbc_x;
            buf[3 * i + 1] = x[j * PAD + 1] + pbc_y;
            buf[3 * i + 2] = x[j * PAD + 2] + pbc_z;
//            fprintf(stderr, "box %d: packing atom %d at %d pbc new pos %e %e %e\n",
//                    this->box_id, j, i, buf[3*i], buf[3*i+1], buf[3*i+2]);
        }
    }
}

void Atom::unpack_comm(int n, int first, double *buf)
{
    int i;

    //#pragma omp for schedule(static)
    for (i = 0; i < n; i++) {
        x[(first + i) * PAD + 0] = buf[3 * i];
        x[(first + i) * PAD + 1] = buf[3 * i + 1];
        x[(first + i) * PAD + 2] = buf[3 * i + 2];
    }

    /* Ensure the new received atoms are not too close to any existing
     * atom in the bin */
    for (i = 0; i < n; i++) {
        check_new_atom(this, first + i, &buf[3 * i]);
    }

    check_ghost_overlap(this);
}

void Atom::pack_reverse(int n, int first, double *buf)
{
    int i;

    //#pragma omp for schedule(static)
    for (i = 0; i < n; i++) {
        buf[3 * i + 0] = f[first + i][0];
        buf[3 * i + 1] = f[first + i][1];
        buf[3 * i + 2] = f[first + i][2];
    }
}

void Atom::unpack_reverse(int n, int *list, double *buf)
{
    int i, j;

    //#pragma omp for schedule(static)
    for (i = 0; i < n; i++) {
        j = list[i];
        f[j][0] += buf[3 * i + 0];
        f[j][1] += buf[3 * i + 1];
        f[j][2] += buf[3 * i + 2];
    }
}

int Atom::pack_border(int i, double *buf, int *pbc_flags)
{
    int m = 0;

    // DSM: Is this branch necessary? pbc_flags are set to 0 if not periodic in a dimension - are the branches not
    // equivalent?
    if (pbc_flags[0] == 0) {
        buf[m++] = x[i * PAD + 0];
        buf[m++] = x[i * PAD + 1];
        buf[m++] = x[i * PAD + 2];
        buf[m++] = type[i]; /* FIXME: Here is a cast from int to double
                               */
    } else {
        buf[m++] = x[i * PAD + 0] + pbc_flags[1] * box.xprd;
        buf[m++] = x[i * PAD + 1] + pbc_flags[2] * box.yprd;
        buf[m++] = x[i * PAD + 2] + pbc_flags[3] * box.zprd;
        buf[m++] = type[i];
    }

    return m;
}

int Atom::unpack_border(int i, double *buf)
{
    if (i == nmax)
        growarray(); // DSM: No change needed for multibox, this grows this box's x,v,f arrays

    int m = 0;
    x[i * PAD + 0] = buf[m++];
    x[i * PAD + 1] = buf[m++];
    x[i * PAD + 2] = buf[m++];
    type[i] = buf[m++];
    return m;
}

int Atom::pack_exchange(int i, double *buf)
{
    int m = 0;
    buf[m++] = x[i * PAD + 0];
    buf[m++] = x[i * PAD + 1];
    buf[m++] = x[i * PAD + 2];
    buf[m++] = v[i * PAD + 0];
    buf[m++] = v[i * PAD + 1];
    buf[m++] = v[i * PAD + 2];
    buf[m++] = type[i];
    return m;
}

int Atom::unpack_exchange(int i, double *buf)
{
    if (i == nmax)
        growarray();

    int m = 0;
    x[i * PAD + 0] = buf[m++];
    x[i * PAD + 1] = buf[m++];
    x[i * PAD + 2] = buf[m++];
    v[i * PAD + 0] = buf[m++];
    v[i * PAD + 1] = buf[m++];
    v[i * PAD + 2] = buf[m++];
    type[i] = buf[m++];
    return m;
}

/* realloc a 2-d double array */

void Atom::sort(Neighbor &neighbor)
{
    fprintf(stderr, "sorting atoms for box %d\n", this->box_id);
    neighbor.binatoms(*this, nlocal);
    //#pragma omp barrier

    binpos = neighbor.bincount;
    bins = neighbor.bins;

    const int atoms_per_bin = neighbor.atoms_per_bin;

    //#pragma omp master
    {
        for (int i = 1; i < neighbor.ntotbins; i++)
            binpos[i] += binpos[i - 1];
        if (copy_size < nmax) {
            destroy_2d_double_array(x_copy);
            destroy_2d_double_array(v_copy);
            destroy_1d_int_array(type_copy);
            x_copy = (double *) create_2d_double_array(nmax, PAD);
            v_copy = (double *) create_2d_double_array(nmax, PAD);
            type_copy = create_1d_int_array(nmax);
            copy_size = nmax;
        }
    }

    //#pragma omp barrier
    double *new_x = x_copy;
    double *new_v = v_copy;
    int *new_type = type_copy;
    double *old_x = x;
    double *old_v = v;
    int *old_type = type;

    //#pragma omp for
    for (int mybin = 0; mybin < neighbor.ntotbins; mybin++) {
        const int start = mybin > 0 ? binpos[mybin - 1] : 0;
        const int count = binpos[mybin] - start;
        for (int k = 0; k < count; k++) {
            const int new_i = start + k;
            const int old_i = bins[mybin * atoms_per_bin + k];

            new_x[new_i * PAD + 0] = old_x[old_i * PAD + 0];
            new_x[new_i * PAD + 1] = old_x[old_i * PAD + 1];
            new_x[new_i * PAD + 2] = old_x[old_i * PAD + 2];
            new_v[new_i * PAD + 0] = old_v[old_i * PAD + 0];
            new_v[new_i * PAD + 1] = old_v[old_i * PAD + 1];
            new_v[new_i * PAD + 2] = old_v[old_i * PAD + 2];
            new_type[new_i] = old_type[old_i];
        }
    }

    //#pragma omp master
    {
        double *x_tmp = x;
        double *v_tmp = v;
        int *type_tmp = type;

        x = x_copy;
        v = v_copy;
        type = type_copy;
        x_copy = x_tmp;
        v_copy = v_tmp;
        type_copy = type_tmp;
    }
    //#pragma omp barrier
}
