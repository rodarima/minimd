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

#include <stdlib.h>
#include "stdio.h"
#include "math.h"
#include "force_lj.h"

ForceLJ::ForceLJ(int ntypes_, int boxes_per_process_)
{
  cutforce = 0.0;
  use_oldcompute = 0;
  reneigh = 1;
  style = FORCELJ;
  ntypes = ntypes_;
  boxes_per_process = boxes_per_process_; // DSM: Multibox change

  cutforcesq = new double[ntypes*ntypes];
  epsilon = new double[ntypes*ntypes];
  sigma6 = new double[ntypes*ntypes];
  sigma = new double[ntypes*ntypes];

  for(int i = 0; i<ntypes*ntypes; i++) {
    cutforcesq[i] = 0.0;
    epsilon[i] = 1.0;
    sigma6[i] = 1.0;
    sigma[i] = 1.0;
  }

  // DSM: Multibox change: one value per box on this process
  eng_vdwl = (double*)malloc(boxes_per_process_ * sizeof(double));    // DSM One of the outputs of compute(). Used in energy()
  virial = (double*)malloc(boxes_per_process_ * sizeof(double));      // DSM One of the outputs of compute(). Used in pressure()
  evflag = (int*)malloc(boxes_per_process_ * sizeof(int));
}
ForceLJ::~ForceLJ() {
  // DSM Multibox
  if (eng_vdwl) { free(eng_vdwl); }
  if (virial) { free(virial); }
  if (evflag) { free(evflag); }
}

void ForceLJ::setup()
{
  for(int i = 0; i<ntypes*ntypes; i++)
    cutforcesq[i] = cutforce * cutforce;
}

// DSM: Decides which compute implementation to call.
// DSM: Simplified function to ensure full neighborlist compute implementation is called regardless of input parameters.
// DSM: "me" is an unused parameter - no MPI calls made in this force implementation
void ForceLJ::compute(Atom &atom, Neighbor &neighbor, Comm &comm, int me)
{
  // DSM: Multibox change. These are now arrays indexed by box id
  eng_vdwl[atom.box_id] = 0;
  virial[atom.box_id] = 0;

  if(evflag[atom.box_id]) {
    return compute_fullneigh(atom, neighbor, me, 1);
  } else {
    return compute_fullneigh(atom, neighbor, me, 0);
  }
}

// DSM: These compute functions are the hotspots according to the 2015 intel webinar:
// https://software.intel.com/en-us/videos/correct-to-correct-and-efficient-a-case-study-with-minimd

//optimised version of compute
//  -MPI + OpenMP (using full neighborlists)
//  -gets rid of fj update (read/write to memory)
//  -use temporary variable for summing up fi
//  -enables vectorization by:
//    -get rid of 2d pointers
//    -use pragma simd to force vectorization of inner loop
//template<int EVFLAG>
void ForceLJ::compute_fullneigh(Atom &atom, Neighbor &neighbor, int me, int EVFLAG)
{
  double t_eng_vdwl = 0;
  double t_virial = 0;

  const int nlocal = atom.nlocal;
  const int nall = atom.nlocal + atom.nghost;
  const double* const x = atom.x;
  double* const f = atom.f;
  const int* const type = atom.type;

  // clear force on own and ghost atoms

  for(int i = 0; i < nlocal; i++) {
    f[i * PAD + 0] = 0.0;
    f[i * PAD + 1] = 0.0;
    f[i * PAD + 2] = 0.0;
  }

  // loop over all neighbors of my atoms
  // store force on atom i

  for(int i = 0; i < nlocal; i++) {
    const int* const neighs = &neighbor.neighbors[i * neighbor.maxneighs];
    const int numneighs = neighbor.numneigh[i];
    const double xtmp = x[i * PAD + 0];
    const double ytmp = x[i * PAD + 1];
    const double ztmp = x[i * PAD + 2];
    const int type_i = type[i];
    double fix = 0;
    double fiy = 0;
    double fiz = 0;

    //pragma simd forces vectorization (ignoring the performance objections of the compiler)
    //also give hint to use certain vectorlength for MIC, Sandy Bridge and WESTMERE this should be be 8 here
    //give hint to compiler that fix, fiy and fiz are used for reduction only

    for(int k = 0; k < numneighs; k++) {
      const int j = neighs[k];
      const double delx = xtmp - x[j * PAD + 0];
      const double dely = ytmp - x[j * PAD + 1];
      const double delz = ztmp - x[j * PAD + 2];
      const int type_j = type[j];
      const double rsq = delx * delx + dely * dely + delz * delz;

      int type_ij = type_i*ntypes+type_j;
      if(rsq < cutforcesq[type_ij]) {
        const double sr2 = 1.0 / rsq;
        const double sr6 = sr2 * sr2 * sr2 * sigma6[type_ij];
        const double force = 48.0 * sr6 * (sr6 - 0.5) * sr2 * epsilon[type_ij];
        fix += delx * force;
        fiy += dely * force;
        fiz += delz * force;

        if(EVFLAG) {
          t_eng_vdwl += sr6 * (sr6 - 1.0) * epsilon[type_ij];
          t_virial += (delx * delx + dely * dely + delz * delz) * force;
        }
      }
    }

    f[i * PAD + 0] += fix;
    f[i * PAD + 1] += fiy;
    f[i * PAD + 2] += fiz;

  }

  t_eng_vdwl *= 4.0;
  t_virial *= 0.5;

  // DSM: Multibox change. These are now arrays
  eng_vdwl[atom.box_id] += t_eng_vdwl;
  virial[atom.box_id] += t_virial;
}

