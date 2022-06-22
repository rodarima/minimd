#ifndef ENDPOINT_H
#define ENDPOINT_H


typedef struct {
    int rank;
    int box;
    int index; /* We only specify how to address the PackBuf */
    int type;
    int dir;
    char name[1024];
} Endpoint;

#include "packbuf.h"
#include "types.h"

int endpoint_is_same(Endpoint *a, Endpoint *b);
void endpoint_get_remote(Sim *sim, Endpoint *local, Endpoint *remote);
void endpoint_set_name(Endpoint *ep);
PackBuf *endpoint_get_packbuf(Sim *sim, Endpoint *ep);

#endif /* ENDPOINT_H */
