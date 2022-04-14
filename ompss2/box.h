#ifndef BOX_H
#define BOX_H

#include "types.h"

void box_realloc(Box *box, int n);
void box_grow_array(Box *box);
void box_add_atom(Box *box, Vec r, Vec v, int type);

#endif /* BOX_H */
