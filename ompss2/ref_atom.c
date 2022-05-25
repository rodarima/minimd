#define _GNU_SOURCE
#define ENABLE_DEBUG 1
#include "types.h"
#include "log.h"
#include "ref.h"
#include "neigh.h"
#include "safe.h"

#include <math.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <errno.h>

typedef struct RefAtom RefAtom;

struct RefAtom {
    int iter;
    int iatom;
    int isghost;
    Vec r;
    Vec f;
    int nnearby;
    RefAtom **nearby;

    /* Matched with the simulation data */
    Box *box;       /* Box where the atom is in the simulation */
    int sim_iatom;  /* Index of the atom in the simulation */
    double dist;    /* Distance to atom (best match) */

    /* Internal flag */
    int nearby_matched;
};

static int
find_position(Box *box, int iref, Vec r, int isghost, Nearby *nearby, double *dist)
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

    dbg("searching for reference atom %d at (%e %e %e) in box %d, iter=%d\n",
            iref, r[X], r[Y], r[Z], box->i, box->iter);

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
        dbg("not found, closest atom %d at (%e %e %e) dist=%e\n",
                jatom, closest[X], closest[Y], closest[Z], minerr);
        return -1;
    } else {
        dbg("found, closest atom %d at (%e %e %e) dist=%e\n",
                jatom, closest[X], closest[Y], closest[Z], minerr);
    }

    return jatom;
}

static void
check_one_nearby(Sim *sim, Box *box, int iatom, Vec refpos)
{
    double dist;
    Nearby *nearby = &box->nearby[iatom];
    int j = find_position(box, -1, refpos, -1, nearby, &dist);

    if (iatom < 0) {
        die("rank%d.box%d: reference nearby atom %d at (%e %e %e) is missing (sim->iter=%d)\n",
                sim->rank, box->i, iatom, refpos[X], refpos[Y], refpos[Z], sim->iter);
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
check_atom(Sim *sim, int iref, Vec refpos, int isghost, int nnearby)
{
    double dist;
    int matches = 0;
    int ibox = -1;
    int iatom;
    Box *box = NULL;

    for (int i = 0; i < sim->nboxes; i++) {
        Box *tmpbox = &sim->box[i];
        int j = find_position(tmpbox, iref, refpos, isghost, NULL, &dist);
        if (j >= 0) {
            if (matches > 0 && !isghost) {
                die("rank%d: iter=%d duplicated atom:\n"
                        "  ref. pos. (%e %e %e)\n"
                        "  box=%d     (%e %e %e)\n"
                        "  box=%d     (%e %e %e)\n",
                        sim->rank, sim->iter,
                        refpos[X], refpos[Y], refpos[Z],
                        ibox,
                        box->r[j][X],
                        box->r[j][Y],
                        box->r[j][Z],
                        i,
                        tmpbox->r[iatom][X],
                        tmpbox->r[iatom][Y],
                        tmpbox->r[iatom][Z]);
            }
            matches++;
            ibox = i;
            iatom = j;
            box = tmpbox;
        }
    }

    if (matches == 0) {
        die("rank%d: reference atom at (%e %e %e) is missing (sim->iter=%d)\n",
                sim->rank, refpos[X], refpos[Y], refpos[Z], sim->iter);
    }

    if (!isghost) {
        Nearby *nearby = &box->nearby[iatom];
        if (nearby->natoms != nnearby) {
            die("reference atom %d at (%e %e %e) mismatch nearby %d (ok=%d)\n",
                    iatom, refpos[X], refpos[Y], refpos[Z],
                    nearby->natoms, nnearby);
        } else {
            dbg("nnearby for atom %d ok\n", iatom);
        }


//        if (ENABLE_REF_NEARBY && !isghost) {
//            check_all_nearby(sim, box, nearby, iatom);
//        }
    }

    dbg("reference atom %d ok, dist=%e\n", iatom, dist);
}

static void
refatom_load_position(Sim *sim, RefAtom **ref, int *n)
{
    int nalloc = 2000; /* Guess */
    RefAtom *refatoms = safe_calloc(nalloc, sizeof(RefAtom));
    int natoms = 0;

    char refpos[PATH_MAX];

    if(snprintf(refpos, PATH_MAX, "%s/atompos.%d.csv",
                sim->refdir, sim->rank) >= PATH_MAX) {
        die("reference atompos path too long\n");
    }

    FILE *reffile;

    if ((reffile = fopen(refpos, "r")) == NULL)
        die("fopen(%s) failed: %s\n", refpos, strerror(errno));

    int refiter;
    int jatom;
    int isghost;
    Vec refr;
    int refnearby;

    /* Skip the header */
    char buf[1024];
    fgets(buf, 1023, reffile);

    int iatom = 0;

    while (1) {
    
        if (natoms == nalloc) {
            /* Grow */
            nalloc *= 2;
            refatoms = safe_realloc(refatoms, nalloc * sizeof(RefAtom));
        }

        RefAtom *a = &refatoms[natoms];

        int ret = fscanf(reffile, "%d,%d,%d,%le,%le,%le,%d\n",
                &a->iter, &a->iatom, &a->isghost,
                &a->r[X], &a->r[Y], &a->r[Z],
                &a->nnearby);

        a->nearby = NULL;
        a->box = NULL;
        a->sim_iatom = -1;
        a->dist = -666.0;

        if (ret == EOF)
            break;
        else if (ret != 7)
            die("wrong number of fields read\n");

        if (a->iter == sim->iter) {
            if (a->iatom != natoms)
                die("inconsistent reference atom index\n");
            /* Only accept it if the iteration and index match */
            natoms++;
        } else if (a->iter > sim->iter) {
            /* Stop as soon as we find an atom of the next iteration */
            break;
        }
    }

    fclose(reffile);

    *ref = refatoms;
    *n = natoms;
}

static void
refatom_load_nearby(Sim *sim, RefAtom *refatoms)
{
    char refpath[PATH_MAX];

    if (snprintf(refpath, PATH_MAX, "%s/nearby.%d.csv",
                sim->refdir, sim->rank) >= PATH_MAX) {
        die("reference nearby path too long\n");
    }

    FILE *reffile;

    if ((reffile = fopen(refpath, "r")) == NULL)
        die("fopen(%s) failed: %s\n", refpath, strerror(errno));

    /* Skip the header */
    char buf[1024];
    fgets(buf, 1023, reffile);

    while (1) {
        int iter, iatom, inearby, jatom;

        int ret = fscanf(reffile, "%d,%d,%d,%d\n",
                &iter, &iatom, &inearby, &jatom);

        if (ret == EOF)
            break;
        else if (ret != 4)
            die("wrong number of fields read\n");

        RefAtom *a = &refatoms[iatom];

        if (a->nearby == NULL) {
            a->nearby = safe_calloc(a->nnearby, sizeof(RefAtom *));
        }

        if (iter < sim->iter)
            continue;

        if (iter > sim->iter)
            break;

        if (inearby >= a->nnearby)
            die("reference atom %d nearby %d exceeds max %d\n",
                    a->iatom, inearby, a->nnearby);

        if (a->nearby[inearby] != NULL)
            die("atom %d attempts to write twice to nearby %d\n",
                    a->iatom, inearby);

        a->nearby[inearby] = &refatoms[jatom];
    }

    fclose(reffile);
}

static void
refatom_load_force(Sim *sim, RefAtom *refatoms)
{
    char refpath[PATH_MAX];

    if (snprintf(refpath, PATH_MAX, "%s/atomforce.%d.csv",
                sim->refdir, sim->rank) >= PATH_MAX) {
        die("reference force path too long\n");
    }

    FILE *reffile;

    if ((reffile = fopen(refpath, "r")) == NULL)
        die("fopen(%s) failed: %s\n", refpath, strerror(errno));

    /* Skip the header */
    char buf[1024];
    fgets(buf, 1023, reffile);

    while (1) {
        int iter, iatom;
        Vec f;

        int ret = fscanf(reffile, "%d,%d,%le,%le,%le\n",
                &iter, &iatom, &f[X], &f[Y], &f[Z]);

        if (ret == EOF)
            break;
        else if (ret != 5)
            die("wrong number of fields read\n");

        if (iter < sim->iter)
            continue;

        if (iter > sim->iter)
            break;

        RefAtom *a = &refatoms[iatom];

        for (int d = X; d <= Z; d++)
            a->f[d] = f[d];
    }

    fclose(reffile);
}

static void
refatom_free(RefAtom *ref, int n)
{
    free(ref);
}

static int
refatom_match(RefAtom *atom, Box *box, double *dist)
{
    int n0, n1;

    if (atom->isghost == 1) {
        n0 = box->nlocal;
        n1 = box->nlocal + box->nghost;
        dbg("searching for reference GHOST atom %d at (%e %e %e) in box %d, iter=%d\n",
                atom->iatom, atom->r[X], atom->r[Y], atom->r[Z], box->i, box->iter);
    } else {
        n0 = 0;
        n1 = box->nlocal;
        dbg("searching for reference LOCAL atom %d at (%e %e %e) in box %d, iter=%d\n",
                atom->iatom, atom->r[X], atom->r[Y], atom->r[Z], box->i, box->iter);
    }

    int found = 0;
    int jatom = -1;
    double curdist;
    double mindist = 1.0;
    Vec closest;
    double maxerr = 1e-6; /* 6 significative figures */
    for (int i = n0; i < n1; i++) {
        Vec *pos = &box->r[i];
        double sqdist = get_distsq(atom->r, *pos);
        curdist = sqrt(sqdist);
        if (curdist < mindist) {
            mindist = curdist;
            jatom = i;
            for (int d = X; d <= Z; d++)
                closest[d] = (*pos)[d];
        }
    }

    double minerr;

    if (dotprod(atom->r) != 0)
        minerr = mindist / sqrt(fabs(dotprod(atom->r)));
    else
        minerr = mindist;

    if (minerr < maxerr) {
        found = 1;
    }

    *dist = minerr;

    if (!found) {
        dbg("ref atom %d NOT FOUND, closest atom %d at (%e %e %e) dist=%e\n",
                atom->iatom, jatom,
                closest[X], closest[Y], closest[Z],
                minerr);

        return -1;
    } else {
        dbg("ref atom %d FOUND, closest atom %d at (%e %e %e) dist=%e\n",
                atom->iatom, jatom,
                closest[X], closest[Y], closest[Z],
                minerr);
    }

    atom->dist = minerr;
    atom->box = box;
    atom->sim_iatom = jatom;

    return jatom;
}

/** Checks one reference atom to match the position computed by the simulation.
 * @param sim       Simulation handle.
 * @param refatom   Base address of reference atoms array.
 * @param irefatom  The index of the atom to be checked.
 * */
static void
refatom_check(Sim *sim, RefAtom *refatom, int irefatom)
{
    double dist, bestdist;
    int matches = 0;
    int ibox = -1;
    int iatom;
    Box *box = NULL;

    RefAtom *a = &refatom[irefatom];

    for (int i = 0; i < sim->nboxes; i++) {
        Box *tmpbox = &sim->box[i];
        int j = refatom_match(a, tmpbox, &dist);
        if (j >= 0) {
            if (matches > 0 && !a->isghost) {
                die("rank%d: iter=%d duplicated atom:\n"
                        "  ref. pos. (%e %e %e)\n"
                        "  box=%d     (%e %e %e)\n"
                        "  box=%d     (%e %e %e)\n",
                        sim->rank, sim->iter,
                        a->r[X], a->r[Y], a->r[Z],
                        ibox,
                        box->r[j][X],
                        box->r[j][Y],
                        box->r[j][Z],
                        i,
                        tmpbox->r[iatom][X],
                        tmpbox->r[iatom][Y],
                        tmpbox->r[iatom][Z]);
            }
            matches++;
            ibox = i;
            iatom = j;
            box = tmpbox;
            bestdist = dist;
        }
    }

    if (matches == 0) {
        die("rank%d: reference atom at (%e %e %e) is missing (sim->iter=%d)\n",
                sim->rank, a->r[X], a->r[Y], a->r[Z], sim->iter);
    }

    if (!a->isghost) {
        Nearby *nearby = &box->nearby[iatom];
        if (nearby->natoms != a->nnearby) {
            die("reference atom %d at (%e %e %e) mismatch nearby %d (ok=%d)\n",
                    iatom, a->r[X], a->r[Y], a->r[Z],
                    nearby->natoms, a->nnearby);
        } else {
            dbg("nnearby for atom %d ok\n", iatom);
        }


//        if (ENABLE_REF_NEARBY && !a->isghost) {
//            check_all_nearby(sim, box, nearby, iatom);
//        }
    }

    dbg("reference atom %d ok, sim atom=%d box=%d dist=%e\n",
            iatom, a->sim_iatom, a->box->i, bestdist);
}

static double
get_reldist(Vec a, Vec b, double maxdist)
{
    if (dotprod(a) != 0 && sqrt(get_distsq(a, b)) > maxdist)
        return sqrt(get_distsq(a, b)) / sqrt(fabs(dotprod(a)));
    else
        return sqrt(get_distsq(a, b));
}

static int
is_close(Vec a, Vec b, double maxrel)
{
    return get_reldist(a, b, maxrel) < maxrel;
}

static void
refatom_nearby_check(Sim *sim, RefAtom *refatoms, int irefatom)
{
    RefAtom *a = &refatoms[irefatom];
    Box *box = a->box;
    Nearby *sim_nearby = &box->nearby[a->sim_iatom];

    for (int i = 0; i < a->nnearby; i++) {
        int jatom = a->nearby[i]->iatom;
        RefAtom *ref_nearby = &refatoms[jatom];

        if (sim_nearby->natoms != a->nnearby) {
            die("sim nearby natoms mismatch\n");
        }

        /* Match the reference nearby atom i by looking at all the reference
         * nearby atoms in the simulation */
        int matches = 0;
        for (int j = 0; j < sim_nearby->natoms; j++) {
            int k = sim_nearby->atom[j];
            if (is_close(box->r[k], ref_nearby->r, 1e-6)) {
                matches++;
            }
        }

        if (matches == 0) {
            die("box%d.atom%d: nearby atom %d not found, pos (%e %e %e)\n",
                    box->i, a->sim_iatom, jatom,
                    a->r[X], a->r[Y], a->r[Z]);
        }

        if (matches > 1) {
            die("box%d.atom%d: nearby atom %d got %d matches\n",
                    box->i, a->sim_iatom, jatom, matches);
        }
    }

    dbg("rank%d: ref atom %d: all %d nearby atoms matched ok\n",
            sim->rank, irefatom, a->nnearby);
}

static void
refatom_force_check(Sim *sim, RefAtom *refatoms, int n)
{

    for (int i = 0; i < n; i++) {
        RefAtom *a = &refatoms[i];

        if (a->isghost)
            continue;

        Box *box = a->box;
        int j = a->sim_iatom;
        Vec f = { box->f[j][X], box->f[j][Y], box->f[j][Z] }; 

        if (!is_close(f, a->f, 1e-6)) {
            die("box%d: reference atom force mismatch:\n"
                "  reference (%e %e %e) i=%d\n"
                "    current (%e %e %e) i=%d\n",
                    box->i,
                    a->f[X], a->f[Y], a->f[Z], a->iatom,
                    f[X], f[Y], f[Z], j);
        } else {
            dbg("box%d: reference atom force ok:\n"
                "  reference (%e %e %e) i=%d\n"
                "    current (%e %e %e) i=%d\n",
                    box->i,
                    a->f[X], a->f[Y], a->f[Z], a->iatom,
                    f[X], f[Y], f[Z], j);
        }
    }

    dbg("rank%d: all forces matched ok\n", sim->rank);
}

static void
compare_atompos(Sim *sim)
{
    RefAtom *refatom;
    int natoms;

    refatom_load_position(sim, &refatom, &natoms);

    refatom_load_nearby(sim, refatom);

    refatom_load_force(sim, refatom);

    /* Check first the ghost atoms */

    for (int i = 0; i < natoms; i++) {
        if (refatom[i].isghost)
            refatom_check(sim, refatom, i);
    }

    for (int i = 0; i < natoms; i++) {
        if (!refatom[i].isghost)
        refatom_check(sim, refatom, i);
    }

//    dbg("--- begin nearby check ---\n");
//
//    for (int i = 0; i < natoms; i++) {
//        if (!refatom[i].isghost)
//            refatom_nearby_check(sim, refatom, i);
//    }

    //refatom_force_check(sim, refatom, natoms);

    refatom_free(refatom, natoms);
}

void
ref_check_atoms(Sim *sim)
{
    if (!ENABLE_REF_ATOMS || sim->refdir == NULL)
        return;

    #pragma oss taskwait /* for debug */

    /* Ensure the iteration numbers are consistent between boxes and the
     * simulation */
    for (int i = 0; i < sim->nboxes; i++) {
        if (sim->iter != sim->box[i].iter) {
            die("rank%d: ref_check_atoms: inconsistent sim->iter=%d != box[%d].iter=%d\n",
                    sim->rank, sim->iter, i, sim->box[i].iter);
        }
    }

    if (ENABLE_REF_ATOMS) {
        compare_atompos(sim);
    }

    dbg("ref_check_atoms ok iter=%d\n", sim->iter);
}
