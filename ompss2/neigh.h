#ifndef NEIGH_H
#define NEIGH_H

#include "types.h"
#include <stdio.h>

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

#endif /* NEIGH_H */
