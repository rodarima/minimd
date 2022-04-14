#include "types.h"


#define BOX_ALLOC_INCR 20000

/* Ensures that at least n atoms fit in the box (both local and ghosts) */
void
box_realloc(Box *box, int n)
{
    if (box->nalloc >= n)
        return;

    box->r = (Vec *) safe_realloc(box->r, n * sizeof(box->r[0]));
    box->v = (Vec *) safe_realloc(box->v, n * sizeof(box->v[0]));
    box->f = (Vec *) safe_realloc(box->f, n * sizeof(box->f[0]));
    box->atomtype = (int *) safe_realloc(box->atomtype, n * sizeof(box->atomtype[0]));
    box->nearby = (Nearby *) safe_realloc(box->nearby, n * sizeof(box->nearby[0]));

    /* Zero the nearby nalloc */
    for (int i = box->nalloc; i < n; i++) {
        Nearby *nearby = &box->nearby[i];
        nearby->nalloc = 0;
        nearby->atom = NULL;
    }

    box->nalloc = n;
}

void
box_grow_array(Box *box)
{
    int n = box->nalloc + BOX_ALLOC_INCR;

    box_realloc(box, n);
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
