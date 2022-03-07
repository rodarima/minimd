#ifndef DOM_H
#define DOM_H

#include "types.h"

/* Returns 1 if the position is inside the domain, otherwise returns 0. */
static int
in_domain(Vec r, Domain dom)
{
    for (int d = X; d <= Z; d++) {
        if (r[d] < dom[d][LO] || r[d] >= dom[d][HI])
            return 0;
    }

    return 1;
}

/* Returns 0 if the position is ouside the domain, and writes the offset
 * in delta. Otherwise returns 1 */
static int
in_domain_delta(Vec r, Domain dom, int delta[NDIM])
{
    int inside = 1;

    for (int d = X; d <= Z; d++) {
        if (r[d] < dom[d][LO]) {
            inside = 0;
            delta[d] = -1;
        } else if (r[d] >= dom[d][HI]) {
            inside = 0;
            delta[d] = +1;
        } else {
            delta[d] = 0;
        }
    }

    return inside;
}

#endif /* DOM_H */
