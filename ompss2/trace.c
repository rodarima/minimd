#define ENABLE_DEBUG 1
#include "log.h"

#include <mpi.h>
#include <stdio.h>
#include <linux/limits.h>

static double hrec = 8.0;
static double hrow = 10.0;
static double tstart = 0.0;
static FILE *tracefile = NULL;
static double xfactor = 10000.0;
static double Nrows = 0;
static int tracing_enabled = 0;

void
trace_open(int rank, int nrows, double duration)
{
    if (!tracing_enabled)
        return;

    double w = duration * xfactor;
    double h = nrows * hrow;
    Nrows = nrows;
    tstart = MPI_Wtime();

    char buf[PATH_MAX];

    if (rank != 0 && rank != 1)
        die("bad rank: %d\n", rank);

    sprintf(buf, "trace.%d.svg", rank);

    tracefile = fopen(buf, "w");

    if (tracefile == NULL)
        die("fopen failed\n");

    fprintf(tracefile, 
    "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
    "<!-- Created with Inkscape (http://www.inkscape.org/) -->\n"
    "<svg\n"
    "   width=\"%dmm\"\n"
    "   height=\"%dmm\"\n"
    "   viewBox=\"0 0 %d %d\"\n"
    "   version=\"1.1\"\n"
    "   id=\"svg5\"\n"
    "   xmlns=\"http://www.w3.org/2000/svg\"\n"
    "   xmlns:svg=\"http://www.w3.org/2000/svg\">\n"
    "  <defs id=\"defs2\" />\n"
    "  <g id=\"layer1\">\n",
        (int) w, (int) h, (int) w, (int) h);



    for (int i = 0; i < nrows + 1; i++) {
        double x0 = 0.0, x1 = w;
        double y0 = i * hrow;
        double y1 = i * hrow;
        fprintf(tracefile,
            "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke=\"#000000\" stroke-width=\"0.2\"/>\n",
            x0, y0, x1, y1);
    }
    fflush(tracefile);
}

void
trace_record_span(double row0, double row1, double t0, double t1, char *name, char *color)
{
    if (!tracing_enabled)
        return;

    if (t1 - t0 < 2e-4)
        t1 = t0 + 2e-4;

    double x = (t0 - tstart) * xfactor;
    double y = row0 * hrow;
    double y1 = row1 * hrow - (hrow - hrec);

    double w = (t1 - t0) * xfactor;
    double h = y1 - y;

    double dur = t1 - t0;

    fprintf(tracefile,
        "<rect fill=\"%s\" \n"
        "    width=\"%f\" height=\"%f\" \n"
        "    x=\"%f\" y=\"%f\" opacity=\"0.6\">\n"
        "  <title>%s\n(dur=%e)</title>\n"
        "</rect>\n", color, w, h, x, y, name, dur);

    fflush(tracefile);
}

void
trace_record(double row, double h, double t0, double t1, char *name, char *color)
{
    if (!tracing_enabled)
        return;

    trace_record_span(row, row + h, t0, t1, name, color);
}

void
trace_event_span(double row0, double row1, char *name, char *color)
{
    if (!tracing_enabled)
        return;

    double t0 = MPI_Wtime();
    double t1 = MPI_Wtime();

    trace_record_span(row0, row1, t0, t1, name, color);
}

void
trace_event(double row, char *name, char *color)
{
    if (!tracing_enabled)
        return;

    trace_event_span(row, row + 1.0, name, color);
}

void
trace_barrier(char *name, char *color)
{
    if (!tracing_enabled)
        return;

    trace_event_span(0.0, (double) Nrows, name, color);
}
