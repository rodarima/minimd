#include "types.h"


#define BOX_ALLOC_INCR 20000

void
box_grow_array(Box *box)
{
    int n = box->nalloc + BOX_ALLOC_INCR;

    box->r = (Vec *) safe_realloc(box->r, n * sizeof(box->r[0]));
    box->v = (Vec *) safe_realloc(box->v, n * sizeof(box->v[0]));
    box->f = (Vec *) safe_realloc(box->f, n * sizeof(box->f[0]));
    box->atomtype = (int *) safe_realloc(box->atomtype, n * sizeof(box->atomtype[0]));
    box->nalloc = n;
}

void
box_add_atom(Box *box, Vec r, Vec v, int type)
{
    int i = box->nlocal;

    if (i == box->nalloc)
        box_grow_array(box);

    for (int d = X; d <= Z; d++) {
        box->r[i][d] = r[d];
        box->v[i][d] = v[d];
    }

    box->atomtype[i] = type;
    box->nlocal++;
}

