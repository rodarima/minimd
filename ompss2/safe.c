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
