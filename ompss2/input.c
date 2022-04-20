/* ----------------------------------------------------------------------
   miniMD is a simple, parallel molecular dynamics (MD) code.   miniMD is
   an MD microapplication in the Mantevo project at Sandia National
   Laboratories ( http://www.mantevo.org ). The primary
   authors of miniMD are Steve Plimpton (sjplimp@sandia.gov) , Paul Crozier
   (pscrozi@sandia.gov) and Christian Trott (crtrott@sandia.gov).

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This library is free software; you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation;
   either version 3 of the License, or (at your option) any later
   version.

   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this software; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA.  See also: http://www.gnu.org/licenses/lgpl.txt .

   For questions, contact Paul S. Crozier (pscrozi@sandia.gov) or
   Christian Trott (crtrott@sandia.gov).

   Please read the accompanying README and LICENSE files.
---------------------------------------------------------------------- */

#include "types.h"

#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#define MAXLINE 4096

static void
safe_fgets(char *s, int size, FILE *stream)
{
    if (fgets(s, size, stream) == NULL) {
        perror("fgets failed");
        abort();
    }
}

static void
usage()
{
    printf(
"-----------------------------------------------------------------------\n"
"%s\n"
"-----------------------------------------------------------------------\n"
"\n"
"miniMD is a simple, parallel molecular dynamics (MD) code,\n"
"which is part of the Mantevo project at Sandia National\n"
"Laboratories ( http://www.mantevo.org ).\n"
"The original authors of miniMD are Steve Plimpton (sjplimp@sandia.gov) ,\n"
"Paul Crozier (pscrozi@sandia.gov) with current\n"
"versions written by Christian Trott (crtrott@sandia.gov).\n"
"\n"
"Commandline Options:\n"
"  -i / --input_file <string>:   set input file to be used\n"
"                                (default: in.lj.miniMD)\n"
"\n"
"  -r / --reference <dir>:       directory with reference outputs\n"
"                                (enable reference comparison)\n"
"\n"
"  -h / --help:                  display this help message\n"
"\n"
"-----------------------------------------------------------------------\n",
        VARIANT_STRING);

    exit(1);
}

void
parse_input_file(Sim *sim, const char *filename)
{
    char line[MAXLINE];
    FILE *fp = fopen(filename, "r");

    if (fp == NULL) {
        fprintf(stderr, "fopen(%s) failed: %s\n", filename,
                strerror(errno));
        abort();
    }

    /* Ignore first two lines (comments) */
    safe_fgets(line, MAXLINE, fp);
    safe_fgets(line, MAXLINE, fp);

    safe_fgets(line, MAXLINE, fp);

    if (strcmp(strtok(line, " \t\n"), "lj") != 0) {
        fprintf(stderr, "Unsupported units '%s', only 'lj' supported\n", line);
        abort();
    }

    safe_fgets(line, MAXLINE, fp);

    if (strcmp(strtok(line, " \t\n"), "none") != 0) {
        fprintf(stderr, "Unsupported datafile '%s', only 'none' supported\n", line);
        abort();
    }

    safe_fgets(line, MAXLINE, fp);

    if (strcmp(strtok(line, " \t\n"), "lj") != 0) {
        fprintf(stderr, "Only 'lj' force type supported\n");
        abort();
    }

    safe_fgets(line, MAXLINE, fp);
    sscanf(line, "%le %le", &sim->epsilon, &sim->sigma);

    safe_fgets(line, MAXLINE, fp);
    sscanf(line, "%d %d %d", &sim->npoints[X], &sim->npoints[Y], &sim->npoints[Z]);

    safe_fgets(line, MAXLINE, fp);
    sscanf(line, "%d", &sim->timesteps);

    safe_fgets(line, MAXLINE, fp);
    sscanf(line, "%le", &sim->dt);

    safe_fgets(line, MAXLINE, fp);
    sscanf(line, "%le", &sim->t_request);

    safe_fgets(line, MAXLINE, fp);
    sscanf(line, "%le", &sim->rho);

    safe_fgets(line, MAXLINE, fp);
    sscanf(line, "%d", &sim->neighbor_period);

    safe_fgets(line, MAXLINE, fp);
    sscanf(line, "%le %le", &sim->R_force, &sim->skin);

    safe_fgets(line, MAXLINE, fp);
    sscanf(line, "%d", &sim->thermo_period);

    safe_fgets(line, MAXLINE, fp);
    sscanf(line, "%d", &sim->nboxes);

    safe_fgets(line, MAXLINE, fp);
    sscanf(line, "%d", &sim->nprocsx);

    safe_fgets(line, MAXLINE, fp);
    sscanf(line, "%d", &sim->nprocsz);

    fclose(fp);

    MPI_Barrier(MPI_COMM_WORLD);
}

static void
check_input(Sim *sim)
{
//    if (sim->nboxes < 3) {
//        fprintf(stderr, "error: at least 3 boxes needed\n");
//        abort();
//    }

    if (sim->timesteps <= 0) {
        fprintf(stderr, "error: timesteps must be > 0\n");
        abort();
    }

    if (sim->npoints[X] == 0 || sim->npoints[Y] == 0 || sim->npoints[Z] == 0) {
        fprintf(stderr, "error: the number of points cannot be 0\n");
        abort();
    }

}

void
parse_input(Sim *sim, int argc, char *argv[])
{
    sim->inputfile = "in.lj.miniMD";
    sim->refdir = NULL;

    /* Skip program name */
    for (int i = 1; i < argc; i++) {
        if ((strcmp(argv[i], "-i") == 0) || (strcmp(argv[i], "--input_file") == 0)) {
            sim->inputfile = argv[++i];
            continue;
        }

        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0)) {
            usage();
            continue;
        }

        if ((strcmp(argv[i], "-r") == 0) || (strcmp(argv[i], "--reference") == 0)) {
            sim->refdir = argv[++i];
            continue;
        }

        fprintf(stderr, "error, unrecognized option '%s'\n", argv[i]);
        usage();
    }

    parse_input_file(sim, (const char *) sim->inputfile);
    check_input(sim);
}
