#include "safe.h"

#include <stdlib.h>
#include <stdio.h>

void *
safe_realloc(void *ptr, size_t size)
{
    void *p = realloc(ptr, size);

    if (p == NULL) {
        perror("realloc failed");
        abort();
    }

    return p;
}

void *
safe_calloc(size_t nmemb, size_t size)
{
    void *p = calloc(nmemb, size);

    if (p == NULL) {
        perror("calloc failed");
        abort();
    }

    return p;
}
