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

#ifndef FORCEEAM_H
#define FORCEEAM_H

#include "atom.hpp"
#include "comm.hpp"
#include "force.hpp"
#include "neighbor.hpp"
#include "threadData.hpp"
#include "types.hpp"
#include <mpi.h>
#include <stdio.h>

class ForceEAM : Force {
public:
    // public variables so USER-ATC package can access them

    double cutmax;

    // potentials as array data

    int nrho, nr;
    int nrho_tot, nr_tot;
    double *frho, *rhor, *z2r;

    // potentials in spline form used for force computation

    double dr, rdr, drho, rdrho;
    double *rhor_spline, *frho_spline, *z2r_spline;

    ForceEAM(int ntypes_);
    virtual ~ForceEAM();
    virtual void compute(Atom &atom, Neighbor &neighbor, Comm &comm, int me);
    virtual void coeff(const char *);
    virtual void setup();
    void init_style();
    double single(int, int, int, int, double, double, double, double &);

    virtual int pack_comm(int n, int iswap, double *buf, int **asendlist);
    virtual void unpack_comm(int n, int first, double *buf);
    int pack_reverse_comm(int, int, double *);
    void unpack_reverse_comm(int, int *, double *);
    double memory_usage();

protected:
    void compute_halfneigh(Atom &atom, Neighbor &neighbor, Comm &comm, int me);
    void compute_fullneigh(Atom &atom, Neighbor &neighbor, Comm &comm, int me);

    // per-atom arrays

    double *rho, *fp;

    int nmax;

    // potentials as file data

    int *map; // which element each atom type maps to

    struct Funcfl {
        char *file;
        int nrho, nr;
        double drho, dr, cut, mass;
        double *frho, *rhor, *zr;
    };
    Funcfl funcfl;

    void array2spline();
    void interpolate(int n, double delta, double *f, double *spline);
    void grab(FILE *, int, double *);

    virtual void read_file(const char *);
    virtual void file2array();

    void bounds(char *str, int nmax, int &nlo, int &nhi);

    void communicate(Atom &atom, Comm &comm);
};

#endif
