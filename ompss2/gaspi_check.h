#ifndef GASPI_CHECK_H
#define GASPI_CHECK_H

#ifndef ENABLE_DEBUG
#define ENABLE_DEBUG 0
#endif

#include "log.h"
#include <GASPI.h>

static void
check_gaspi(gaspi_return_t ret, const char *fn, const char *file, int line)
{
    if (ret == GASPI_SUCCESS)
        return;

    char *msg;
    gaspi_print_error(ret, &msg);
    die("%s (%s at %s:%d)\n", msg, fn, file, line);
}

#define CHECK(x) check_gaspi(x, #x, __FILE__, __LINE__)

#endif /* GASPI_CHECK_H */
