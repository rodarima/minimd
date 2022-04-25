#include "types.h"
#include "neigh.h"
#include "dom.h"

#define BIN_ALLOC_INCR 100
#define NEARBY_ALLOC_INCR 100

static void
clear_bins(Box *box)
{
    for (int i = 0; i < box->nbinsalloc; i++) {
        Bin *bin = &box->bin[i];
        bin->natoms = 0;
    }
}

static void
add_atom_bin(Bin *bin, int iatom)
{
    if (bin->natoms == bin->nalloc) {
        int n = bin->nalloc + BIN_ALLOC_INCR * sizeof(int);
        bin->atom = (int *) safe_realloc(bin->atom, (size_t) n);
        bin->nalloc = n;
    }
    bin->atom[bin->natoms++] = iatom;
}

static void
add_nearby_atom(Nearby *nearby, int iatom)
{
    if (nearby->natoms == nearby->nalloc) {
        int n = nearby->nalloc + NEARBY_ALLOC_INCR * sizeof(int);
        nearby->atom = (int *) safe_realloc(nearby->atom, (size_t) n);
        nearby->nalloc = n;
    }
    nearby->atom[nearby->natoms++] = iatom;
}

#pragma oss task label("build_nearby_atoms_box") \
    in(box->r) \
    in(box->bin) \
    out(box->nearby)
static void
build_nearby_atoms_box(Sim *sim, Box *box)
{
    dbg("building nearby atoms for box %d\n", box->i);
    /* Build nearby lists only for local atoms */
    for (int iatom = 0; iatom < box->nlocal; iatom++) {
        Vec ri = { box->r[iatom][X], box->r[iatom][Y], box->r[iatom][Z] };
        int itype = box->atomtype[iatom];

        /* Clear the previous nearby list */
        Nearby *nearby = &box->nearby[iatom];
        nearby->natoms = 0;

        /* Find the current atom bin */
        int iindbin = get_atom_bin(sim, box, box->r[iatom]);
        Bin *ibin = &box->bin[iindbin];

        /* Find nearing atoms in the stencil bins */
        for (int j = 0; j < box->nstencil; j++) {
            int jindbin = iindbin + box->stencil[j];

            /* The bin should not fall outside the bin range, even if we
             * iterate through atoms in the corner, as the local atoms
             * are inside the box domain and the stencil should only
             * reach bins within the halo domain. */
            if (jindbin < 0 || jindbin >= box->nbinsalloc) {
                /* FIXME: we shouldn't enter here */
                continue;
                fprintf(stderr, "near atom bin %d is outside the range\n", jindbin);
                abort();
            }

            Bin *jbin = &box->bin[jindbin];

            /* Find atoms within R_neigh in the bin */
            for (int k = 0; k < jbin->natoms; k++) {
                int jatom = jbin->atom[k];

                /* Skip if we are find the current atom */
                if (iatom == jatom) {
                    if (ibin != jbin) {
                        fprintf(stderr, "same atom found in different bins\n");
                        abort();
                    }
                    continue;
                }

                Vec rj = { box->r[jatom][X], box->r[jatom][Y], box->r[jatom][Z] };
                int jtype = box->atomtype[jatom];

                double dist_sq = get_distsq(ri, rj);
                int pairtype = itype * sim->ntypes + jtype;
                double R_neigh_sq = sim->R_neigh_sq[pairtype];

                /* You have gone too far */
                if (dist_sq >= R_neigh_sq)
                    continue;

                add_nearby_atom(nearby, jatom);
            }
        }
    }

    box->fresh_ghost = 1;
}

#pragma oss task label("bin_atoms") \
    in(box->r) \
    out(box->bin)
static void
bin_atoms(Sim *sim, Box *box)
{
    clear_bins(box);

    for (int j = 0; j < box->nlocal + box->nghost; j++) {
        int ibin = get_atom_bin(sim, box, box->r[j]);
        Bin *bin = &box->bin[ibin];

        add_atom_bin(bin, j);
    }
}

/* Builds the list of nearby atoms for each local atom of the boxes */
void
build_nearby_atoms(Sim *sim)
{
    for (int i = 0; i < sim->nboxes; i++) {
        Box *box = &sim->box[i];
        bin_atoms(sim, box);
        build_nearby_atoms_box(sim, box);
    }
}
