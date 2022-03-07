#include "types.h"
#include "neigh.h"
#include "dom.h"

static int
get_bin_index(Box *box, int i[NDIM])
{
    int index = i[Z] * box->nbinshalo[Y] * box->nbinshalo[X]
        + i[Y] * box->nbinshalo[X]
        + i[X];

    if (index < 0 || index >= box->nbinsalloc)
        abort();

    return index;
}

/* convert xyz atom coords into local bin #
   take special care to insure ghost atoms with
   coord >= prd or coord < 0.0 are put in correct bins */

static inline int
get_atom_bin(Sim *sim, Box *box, Vec r)
{
    int i[NDIM];

    /* The atom position must be inside the halo domain */
    if (!in_domain(r, box->domhalo)) {
        fprintf(stderr, "atom outside halo domain: %e %e %e\n",
                r[X], r[Y], r[Z]);
        abort();
    }

    for (int d = X; d <= Z; d++) {
        /* Compute the relative position of the atom inside the halo
         * domain, and then just obtain the bin index, dividing by the
         * bin length */
        double delta = r[d] - box->domhalo[d][LO];

        /* FIXME: multiply by the inverse to avoid expensive division */
        i[d] = delta / sim->binlen[d];
    }

    return get_bin_index(box, i);
}

static void
bin_atoms(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        for (int j = 0; j < box->nlocal + box->nghost; j++) {
            int bin = get_atom_bin(sim, box, box->r[j]);
            //fprintf(stderr, "atom %d at %e %e %e gets bin %d\n",
            //        j, box->r[j][X], box->r[j][Y], box->r[j][Z], bin);
        }
    }

}

/* Builds the list of neighbor atoms in each local atom of the boxes */
void
build_neighlist(Sim *sim)
{
    bin_atoms(sim);
}
