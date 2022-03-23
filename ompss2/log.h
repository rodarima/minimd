#ifndef LOG_H
#define LOG_H

#include <stdio.h>
#include <stdlib.h>

 /* Debug macros */
#ifdef ENABLE_DEBUG
# define dbg(...) fprintf(stderr, __VA_ARGS__);
#else
# define dbg(...)
#endif

#define err(...) fprintf(stderr, __VA_ARGS__);
#define die(...) do { err("fatal: " __VA_ARGS__); abort(); } while (0) 

#endif /* LOG_H */
