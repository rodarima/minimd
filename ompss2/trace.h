#ifndef TRACE_H
#define TRACE_H

void
trace_open(int rank, int nrows, double duration);

void
trace_record(double row, double h, double t0, double t1, char *name, char *color);

void
trace_record_span(double row0, double row1, double t0, double t1, char *name, char *color);

void
trace_event(double row, char *name, char *color);

void
trace_event_span(double row0, double row1, char *name, char *color);

void
trace_barrier(char *name, char *color);

#endif /* TRACE_H */
