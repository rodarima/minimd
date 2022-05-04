#ifndef SAFE_H
#define SAFE_H

#include <stdlib.h>

void *safe_realloc(void *ptr, size_t size);
void *safe_calloc(size_t nmemb, size_t size);

#endif /* SAFE_H */
