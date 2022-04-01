#ifndef NEIGH_H
#define NEIGH_H

#define ENABLE_DEBUG 0
#include "types.h"
#include "log.h"
#include "dom.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static inline int
opposite_neigh(int i)
{
	return (NNEIGH - 1) - i;
}

static inline int
delta2neigh(int delta[NDIM])
{
    if (delta[X] == 0 && delta[Y] == 0 && delta[Z] == 0)
        return -1;

    int nd = NNEIGHDIM;
    int ns = NNEIGHSIDE;

    /* Transform delta into zero based offsets */
    int off[NDIM] = {
        delta[X] + ns,
        delta[Y] + ns,
        delta[Z] + ns
    };

    int index = off[Z] * nd * nd + off[Y] * nd + off[X];

    /* Skip the center box at (0,0,0) */
    if (index >= NNEIGH / 2)
        index--;

    //fprintf(stderr, "delta2neighneigh: delta(%2d %2d %2d) -> %d\n", 
    //	    delta[X], delta[Y], delta[Z], index);

    return index;
}

static inline int
delta2subdom(int delta[NDIM])
{
    /* Transform delta into zero based offsets */
    int off[NDIM] = {
        delta[X] + 1,
        delta[Y] + 1,
        delta[Z] + 1
    };

    return off[Z] * 3 * 3 + off[Y] * 3 + off[X];
}

static int
get_bin_index(Box *box, int i[NDIM])
{
    int index = i[Z] * box->nbinshalo[Y] * box->nbinshalo[X]
        + i[Y] * box->nbinshalo[X]
        + i[X];

    if (index < 0 || index >= box->nbinsalloc)
        abort();

    return index;
}

/* convert xyz atom coords into local bin #
   take special care to insure ghost atoms with
   coord >= prd or coord < 0.0 are put in correct bins */

static inline int
get_atom_bin(Sim *sim, Box *box, Vec r)
{
    int i[NDIM];

    /* The atom position must be inside the halo domain */
    if (ENABLE_DOMAIN_CHECK && !in_domain(r, box->domhalo)) {
        die("atom outside halo domain: %e %e %e\n",
                r[X], r[Y], r[Z]);
    }

    for (int d = X; d <= Z; d++) {
        /* Compute the relative position of the atom inside the halo
         * domain, and then just obtain the bin index, dividing by the
         * bin length */
        double delta = r[d] - box->domhalo[d][LO];
        i[d] = delta * sim->invbinlen[d];
    }

    return get_bin_index(box, i);
}

static double
dotprod(Vec v)
{
    double sum = 0.0;

    for (int d = X; d <= Z; d++) {
        sum += v[d] * v[d];
    }

    return sum;
}

static void
vec_diff(Vec delta, Vec ri, Vec rj)
{
    for (int d = X; d <= Z; d++)
        delta[d] = ri[d] - rj[d];
}

static double
get_distsq(Vec ri, Vec rj)
{
    Vec delta;
    vec_diff(delta, ri, rj);
    return dotprod(delta);
}

#endif /* NEIGH_H */
