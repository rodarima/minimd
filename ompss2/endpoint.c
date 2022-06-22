#define ENABLE_DEBUG 0
#include "log.h"
#include "endpoint.h"
#include "packbuf.h"

#include <stddef.h>
#include <stdio.h>
#include <string.h>

enum pb_dir
reverse_dir(enum pb_dir dir)
{
    if (dir == PB_SEND)
        return PB_RECV;
    else
        return PB_SEND;
}

//static void
//get_box_coords(Sim *sim, Endpoint *ep, int *coords)
//{
//    int rank = ep->rank;
//    int ibox = ep->box;
//
//    /* We need a global table of coordinates for all boxes */
//}
//
///* This is the only function that properly computes the opposite endpoint. It
// * may happen that both endpoints are not on the current rank. */
//void
//endpoint_get_other_side(Sim *sim, Endpoint *local, Endpoint *remote)
//{
//    Box *box = &sim->box[local->box];
//    Neigh *neigh = &box->neigh[local->index];
//
//    remote->rank = neigh->rank;
//    remote->box = neigh->boxid;
//    remote->index = local->index;
//    remote->type = local->type;
//    remote->dir = reverse_dir(local->dir);
//}

static void
get_remote_destination(Sim *sim, Endpoint *src, Endpoint *dst)
{
    if (src->dir != PB_SEND)
        die("wrong direction for source\n");

    /* The send direction is very easy, just follow the neighbor box
     * at index to know which is the destination box */
    Box *srcbox = &sim->box[src->box];
    Neigh *neigh = &srcbox->neigh[src->index];

    dst->rank = neigh->rank;
    dst->box = neigh->boxid;
    dst->index = src->index;
    dst->type = src->type;
    dst->dir = reverse_dir(src->dir);
}

static void
get_remote_source(Sim *sim, Endpoint *dst, Endpoint *src)
{
    if (dst->dir != PB_RECV)
        die("wrong direction for destination\n");

    /* To get the source, we invert the neighbor */
    Box *dstbox = &sim->box[dst->box];
    Neigh *dstneigh = dstbox->neigh[src->index].opposite;

    src->rank = dstneigh->rank;
    src->box = dstneigh->boxid;
    src->index = dst->index;
    src->type = dst->type;
    src->dir = reverse_dir(dst->dir);
}

/** Computes the remote side of a *local* endpoint. The local endpoint must
 * reside in the current rank. */
void
endpoint_get_remote(Sim *sim, Endpoint *local, Endpoint *remote)
{
    if (local->rank != sim->rank)
        die("the local endpoint doesn't belong to this rank\n");

    if (local->dir == PB_SEND)
        get_remote_destination(sim, local, remote);
    else
        get_remote_source(sim, local, remote);
}

PackBuf *
endpoint_get_packbuf(Sim *sim, Endpoint *ep)
{
    if (ep->rank != sim->rank)
        return NULL;

    if (ep->box < 0 || ep->box >= sim->nboxes)
        return NULL;

    if (ep->index < 0 || ep->index >= NNEIGH)
        return NULL;

    if (ep->type < 0 || ep->type >= PB_NTYPES)
        return NULL;

    if (ep->dir < 0 || ep->dir >= PB_NDIR)
        return NULL;

    Box *box = &sim->box[ep->box];
    Neigh *neigh = &box->neigh[ep->index];
    PackBuf *pb = &neigh->pb[ep->type][ep->dir];

    return pb;
}

int
endpoint_is_same(Endpoint *a, Endpoint *b)
{

    /* Skip the name part of the endpoint, as it may contain garbage
     * after the string ends */
    if (memcmp(a, b, offsetof(Endpoint, name)) == 0)
        return 1;

    return 0;
}

void
endpoint_set_name(Endpoint *ep)
{
    sprintf(ep->name, "Endpoint(rank=%d box=%d index=%d type=%s dir=%s)",
            ep->rank, ep->box, ep->index,
            PB_TYPENAME(ep->type),
            PB_DIRNAME(ep->dir));
}
