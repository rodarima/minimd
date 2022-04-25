#define ENABLE_DEBUG 1

#include "types.h"
#include "log.h"
#include "neigh.h"

#include <linux/limits.h>
#include <errno.h>
#include <string.h>
#include <math.h>

static int
find_position(Box *box, Vec r, int isghost, Nearby *nearby, double *dist)
{
    int n0, n1;
    int *list = NULL;

    if (isghost == 1) {
        n0 = box->nlocal;
        n1 = box->nlocal + box->nghost;
    } else if (isghost == 0) {
        n0 = 0;
        n1 = box->nlocal;
    } else {
        n0 = 0;
        n1 = nearby->natoms;
        list = nearby->atom;
    }

    int found = 0;
    int jatom = -1;
    double curdist;
    double mindist = 1.0;
    Vec closest;
    double maxerr = 1e-6; /* 6 significative figures */
    for (int i = n0; i < n1; i++) {
        int j;
        if (list) {
            j = list[i];
        } else {
            j = i;
        }
        Vec *pos = &box->r[j];
        double sqdist = get_distsq(r, *pos);
        curdist = sqrt(sqdist);
        if (curdist < mindist) {
            mindist = curdist;
            jatom = i;
            for (int d = X; d <= Z; d++)
                closest[d] = (*pos)[d];
        }
    }

    double minerr;

    if (dotprod(r) != 0)
        minerr = mindist / sqrt(fabs(dotprod(r)));
    else
        minerr = mindist;

    if (minerr < maxerr) {
        found = 1;
    }

    *dist = minerr;

    if (!found) {
        dbg("not found, closest atom %d at (%e %e %e) [to (%e %e %e)]\n",
                jatom, closest[X], closest[Y], closest[Z],
                r[X], r[Y], r[Z]);
        return -1;
    }

    return jatom;
}

static void
check_one_nearby(Sim *sim, Box *box, int iatom, Vec refpos)
{
    double dist;
    Nearby *nearby = &box->nearby[iatom];
    int j = find_position(box, refpos, -1, nearby, &dist);

    if (iatom < 0) {
        die("reference nearby atom %d at (%e %e %e) is missing (sim->iter=%d)\n",
                iatom, refpos[X], refpos[Y], refpos[Z], sim->iter);
    }

    //dbg("reference nearby %d for atom %d ok, dist=%e\n", j, iatom, dist);
}

static void
check_all_nearby(Sim *sim, Box *box, Nearby *nearby, int iatom)
{
    char path[PATH_MAX];

    if(snprintf(path, PATH_MAX, "%s/atomneigh.csv", sim->refdir) >= PATH_MAX) {
        die("reference atomneigh.csv path too long\n");
    }

    FILE *ref;

    if ((ref = fopen(path, "r")) == NULL)
        die("fopen(%s) failed: %s\n", path, strerror(errno));

    int refiter;
    int jatom, jneigh, jatomneigh;
    Vec refr;

    /* Skip the header */
    char buf[1024];
    fgets(buf, 1023, ref);

    /* Iterate through all reference nearby atoms and ensure they exist in the
     * simulation (instead of the other way around). */
    while(1) {
        int ret = fscanf(ref, "%d,%d,%d,%d,%le,%le,%le\n",
                &refiter, &jatom, &jneigh, &jatomneigh,
                &refr[X], &refr[Y], &refr[Z]);

        if (ret == EOF)
            break;
        else if (ret != 7)
            die("wrong number of fields read\n");

        if (refiter == box->iter && iatom == jatom) {
            check_one_nearby(sim, box, iatom, refr);
        } else if (refiter > box->iter) {
            break;
        }
    }

    fclose(ref);

}

static void
check_atom(Sim *sim, Box *box, int iatomref, Vec refpos,
        int isghost, int nnearby)
{
    double dist;

    int iatom = find_position(box, refpos, isghost, NULL, &dist);

    if (iatom < 0) {
        die("reference atom %d at (%e %e %e) is missing (sim->iter=%d)\n",
                iatom, refpos[X], refpos[Y], refpos[Z], sim->iter);
    }

    Nearby *nearby = &box->nearby[iatom];
    if (nearby->natoms != nnearby) {
        die("reference atom %d at (%e %e %e) mismatch nearby %d (ok=%d)\n",
                iatom, refpos[X], refpos[Y], refpos[Z],
                nearby->natoms, nnearby);
    }

    dbg("reference atom %d ok, dist=%e\n", iatom, dist);

    if (ENABLE_REF_NEARBY && !isghost) {
        check_all_nearby(sim, box, nearby, iatom);
    }
}

static void
compare_atompos(Sim *sim)
{
    if (sim->nboxes > 1)
        die("only 1 box supported in reference mode\n");

    if (sim->nranks > 1)
        die("only 1 rank supported in reference mode\n");

    Box *box = &sim->box[0];
    char refpos[PATH_MAX];

    if(snprintf(refpos, PATH_MAX, "%s/atompos.csv",
                sim->refdir) >= PATH_MAX) {
        die("reference atompos.csv path too long\n");
    }

    FILE *ref;

    if ((ref = fopen(refpos, "r")) == NULL)
        die("fopen(%s) failed: %s\n", refpos, strerror(errno));

    int refiter;
    int jatom;
    int isghost;
    Vec refr;
    int refnearby;

    /* Skip the header */
    char buf[1024];
    fgets(buf, 1023, ref);

    int iatom = 0;

    while(1) {
        int ret = fscanf(ref, "%d,%d,%d,%le,%le,%le,%d\n",
                &refiter, &jatom, &isghost,
                &refr[X], &refr[Y], &refr[Z],
                &refnearby);

        if (ret == EOF)
            break;
        else if (ret != 7)
            die("wrong number of fields read\n");

        if (refiter == box->iter) {
            check_atom(sim, box, iatom, refr, isghost, refnearby);
            iatom++;
        } else if (refiter > box->iter) {
            break;
        }
    }

    fclose(ref);
}

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
    in(sim->Ekin, sim->Epot, sim->Etot)
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

    while(1) {
        int ret = fscanf(ref, "%d,%lf,%lf,%lf\n",
                &refiter, &refEpot, &refEkin, &refEtot);

        if (ret == EOF)
            break;
        else if (ret != 4)
            die("wrong number of fields read\n");

        if (refiter != sim->iter - 1)
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
        dbg("ref energy ok\n");
}

void
ref_check_atoms(Sim *sim)
{
    if (!ENABLE_REF_ATOMS || sim->refdir == NULL)
        return;

    #pragma oss taskwait /* for debug */

    if (ENABLE_REF_ATOMS) {
        compare_atompos(sim);
    }

    dbg("ref comp ok\n");
}

void
ref_check_energy(Sim *sim)
{
    if (!ENABLE_REF_ENERGY || sim->refdir == NULL)
        return;

    if (ENABLE_REF_ENERGY) {
        check_energy(sim);
    }
}
