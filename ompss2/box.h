#ifndef BOX_H
#define BOX_H

#include "types.h"

void box_realloc(Box *box, int n);
void box_grow_array(Box *box);
void box_add_atom(Box *box, Vec r, Vec v, int type);

#include "packbuf.h"

void box_packbuf_clear(Box *box, enum pb_type type, enum pb_dir dir);
void box_packbuf_switch(Box *box, enum pb_type type, enum pb_dir dir, int from, int to);

#endif /* BOX_H */
