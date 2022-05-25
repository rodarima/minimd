#define ENABLE_DEBUG 0

#include "types.h"
#include "log.h"
#include "neigh.h"

#include <linux/limits.h>
#include <errno.h>
#include <string.h>
#include <math.h>


static int
too_far(double val, double ref, double maxrel)
{
    double abserr = fabs(val - ref);
    double relerr = abserr / fabs(ref);

    if (relerr > maxrel)
        return 1;
    
    return 0;
}

#pragma oss task label("check_energy") \
    in(sim->Ekin, sim->Epot, sim->Etot) \
    in(sim->iter)
static void
check_energy(Sim *sim)
{
    if (sim->rank != 0)
        return;

    char path[PATH_MAX];

    if(snprintf(path, PATH_MAX, "%s/energy.csv", sim->refdir) >= PATH_MAX) {
        die("reference energy.csv path too long\n");
    }

    FILE *ref;

    if ((ref = fopen(path, "r")) == NULL)
        die("fopen(%s) failed: %s\n", path, strerror(errno));

    /* Skip the header */
    char buf[1024];
    fgets(buf, 1023, ref);

    double maxrel = 1e-6;
    int tested = 0;

    int refiter;
    double refEkin, refEpot, refEtot;

    while (1) {
        int ret = fscanf(ref, "%d,%lf,%lf,%lf\n",
                &refiter, &refEpot, &refEkin, &refEtot);

        if (ret == EOF)
            break;
        else if (ret != 4)
            die("wrong number of fields read\n");

        if (refiter != sim->iter)
            continue;

        tested = 1;

        if (too_far(sim->Ekin, refEkin, maxrel))
            die("iter %d: kinetic energy mismatch, value=%e reference=%e\n",
                    sim->iter, sim->Ekin, refEkin);

        if (too_far(sim->Epot, refEpot, maxrel))
            die("iter %d: potential energy mismatch, value=%e reference=%e\n",
                    sim->iter, sim->Epot, refEpot);

        if (too_far(sim->Etot, refEtot, maxrel))
            die("iter %d: total energy mismatch, value=%e reference=%e\n",
                    sim->iter, sim->Etot, refEtot);

        break;
    }

    fclose(ref);

    if (!tested)
        die("energy for iteration %d not found in reference file '%s'\n",
                sim->iter, path);
    else
        dbg("ref_check_energy ok iter=%d\n", sim->iter);
}

void
ref_check_energy(Sim *sim)
{
    if (ENABLE_REF_ENERGY) {
        if (sim->refdir == NULL)
            return;

        check_energy(sim);
    }
}
