#ifndef HIST_H
#define HIST_H

#include <stdlib.h>
#include <stdio.h>

static void
hist_clear(Hist *h)
{
    for (int i = 0; i < h->nbins; i++) {
        h->count[i] = 0;
    }
}

static void
hist_init(Hist *h, int nbins, const char *filepath, double delta)
{
    if (nbins > HIST_MAX_NBINS)
        abort();

    h->filepath = filepath;
    h->delta = delta;
    h->nbins = nbins;

    hist_clear(h);

    /* Truncate the file and print the header */
    FILE *f = fopen(filepath, "w");
    fprintf(f, "iter,bin,limit,count\n");
    fclose(f);
}

static void
hist_add(Hist *h, double value)
{
    int bin = value / h->delta;

    if (bin >= h->nbins)
        bin = h->nbins - 1;
    else if (bin < 0)
        abort();

    h->count[bin]++;
}

static void
hist_print(Hist *h, int iter)
{
    FILE *f = fopen(h->filepath, "a");
    for (int i = 0; i < h->nbins; i++) {
        double limit = i * h->delta;
        fprintf(f, "%d,%d,%e,%d\n", iter, i, limit, h->count[i]);
    }
    fclose(f);
}

#endif /* HIST_H */
