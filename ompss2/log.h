#ifndef LOG_H
#define LOG_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#ifndef ENABLE_DEBUG
#error ENABLE_DEBUG must be defined
#endif

#define ENABLE_GLOBAL_DEBUG 1

 /* Debug macros */
# define dbg(...) do { \
    if(ENABLE_DEBUG && ENABLE_GLOBAL_DEBUG) fprintf(stderr, __VA_ARGS__); \
} while (0)

#define err(...) \
    fprintf(stderr, __VA_ARGS__);

#define die(...) do { \
    err("fatal: " __VA_ARGS__); \
    fflush(stderr); \
    if (ENABLE_SLOW_DEATH) { \
        for (int __die__i = 0; __die__i < DEATH_SLEEP/3; __die__i++) { \
            sleep(3); \
            err("fatal (repeat): " __VA_ARGS__); \
        } \
    } \
    abort(); \
} while (0) 

#endif /* LOG_H */
