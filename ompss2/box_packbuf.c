#include "box.h"

void
box_packbuf_clear(Box *box, enum pb_type type, enum pb_dir dir)
{
    for (int i = 0; i < NNEIGH; i++) {
        PackBuf *pb = box->pb[type][dir][i];
        packbuf_clear(pb);
    }
}

void
box_packbuf_switch(Box *box, enum pb_type type, enum pb_dir dir,
        int from, int to)
{
    for (int i = 0; i < NNEIGH; i++) {
        PackBuf *pb = box->pb[type][dir][i];
        packbuf_debug_switch(pb, from, to);
    }
}
