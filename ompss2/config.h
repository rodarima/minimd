#ifndef CONFIG_H
#define CONFIG_H

/* Debug checks */

#define MIN_DISTSQ 0.03

/* Extra checks: these allow early catch of problems but introduce large
 * or very large overheads, so they should be disabled when taking into
 * account the performance of the simulation. */

/* Print a histogram of the force magnitudes per box. It should be
 * smooth. */
#define ENABLE_FHIST 0

/* Print a histogram of the velocity magnitudes per box. */
#define ENABLE_VHIST 0

/* Print a histogram of the distance between nearby atoms. */
#define ENABLE_DHIST 0

/* Compute the energy during the simulation. Needed to validate the
 * results. */
#define ENABLE_REALTIME_ENERGY 1

/* Halts the simulation if the force is too large */
#define ENABLE_MAX_FORCE_CHECK 0
#define MAX_FORCE_SQ (5000*5000)

/* Ensures the atom doesn't move more than the size of the box in a
 * single iteration */
#define ENABLE_MAX_VELOCITY_CHECK 0

/* Writes the position of the atoms per iteration. Introduces a lot of
 * overhead */
#define ENABLE_ATOM_TRACKING 0

/* Halts the simulation if an atom doesn't interact with at least half
 * the neighbors (they are too far away to interact). This may happen
 * with too many time steps without re-neighboring. */
#define ENABLE_MIN_INTERACTIONS_CHECK 1

/* Checks that nearby atoms don't move more than a set limit. */
#define ENABLE_MAX_JUMP_CHECK 1

/* Counts the number of total force interactions */
#define ENABLE_COUNT_INTERACTIONS 0

/* Enable domain checks: ensures the atoms are inside the box or other
 * space domains */
#define ENABLE_DOMAIN_CHECK 0

/* Checks the number of atoms is expected before and after an operation.
 * Needs task wait so it can cause other bugs to disappear. */
#define ENABLE_ATOM_COUNT_CHECK 0

/* Uses only these many atoms. Use 0 to run normally */
#define ENABLE_ONLY_NTOTATOMS 0

/* Ensure that no new atom is too close to a local atom (slow) */
//#define ENABLE_NEW_ATOM_CHECK

/* Ensure that no ghost atom is too close to a local atom (slow) */
//#define ENABLE_GHOST_ATOM_CHECK

/* Performs a ping pong test to check communications (slow) */
#define ENABLE_COMM_TEST 0

/* Update the force following bin order instead of atom sequence */
#define ENABLE_FORCE_BY_BINS 0

/* Use MPI_Isend and MPI_Irecv */
#define ENABLE_NONBLOCKING_MPI 1

/* Use the non-blocking mode of TAMPI, which blocks the release of the task
 * until the MPI requests have been completed */
#define ENABLE_NONBLOCKING_TAMPI 0

/* Use MPI_Waitall if needed */
#define NEED_EXPLICIT_WAIT (ENABLE_NONBLOCKING_MPI && !ENABLE_NONBLOCKING_TAMPI)

/* Use MPI_Wait instead of MPI_Waitall to individually wait for each
 * MPI request (useful for debugging deadlocks). */
#define ENABLE_SEQUENTIAL_MPIWAIT 1

/* If enabled, tagaspi will be used to exchange ghost positions */
#define ENABLE_GASPI 0

/* Wait a large delay before aborting when a problem occurs, so a
 * debugger can be attached. Also allows other aborts to trip. */
#define ENABLE_SLOW_DEATH 0

/* Sleep for a day before aborting (if enabled) */
#define DEATH_SLEEP (3600*24)

/* Track the state of the PackBuf to detect concurrent access */
#define ENABLE_PACKBUF_STATE 1

/* Same but for external usage */
#define ENABLE_PACKBUF_DEBUG_STATE 1

/* Compare energy with reference (needs energy.csv). */
#define ENABLE_REF_ENERGY 1 /* Default 1 */

/* Compare position of all atoms with reference (needs atompos.csv).
 * Warning: slow. */
#define ENABLE_REF_ATOMS 0 /* Default 0 */

/* Compare all nearby atoms with reference (needs atomneigh.csv).
 * Warning: very slow. */
#define ENABLE_REF_NEARBY 0 /* Default 0 */

/* -------------------- DANGER ZONE BEGINS -------------------------- */

/* These options cause the energy values reported by the simulation to
 * diverge from the reference version, thus preventing a direct
 * comparison of the values. This corrections are followed from the
 * book Understanding molecular simulation: from algorithms to
 * applications, Daan Frenkel and Berend Smit (2002, Academic Press).
 * They cause the total energy to be closer to 0 with no drift, so
 * better error checks can be done. Disable them if you need to compare
 * the energy values with the reference. These don't affect the
 * position, velocity or forces in the atoms during the simulation. */

/* Correct the potential energy at R_force for atoms that leave the
 * interaction zone (also referred to e_cut) */
#define ENABLE_ECUT_CORRECTION 0

/* Correct the kinetic energy scale by using ntotatoms instead of
 * ntotatoms - 1, as was being done in the reference. Leads to much
 * smaller errors in total energy when active, but breaks compatibility
 * with reference version values. */
#define ENABLE_NTOTATOMS_CORRECTION 0

/* When the two corrections are enabled, we can measure the total energy
 * at the end of the simulation and check if it diverges. Without the
 * corrections it drifts, so the check will always fail. */
#define MAX_ENERGY_REL_ERROR 0.0005

/* -------------------- END OF DANGER ZONE -------------------------- */

/* ----------------------- Check pre-requisites -------------------- */

#if ENABLE_REF_COMPARE
# if !ENABLE_ATOM_TRACKING
#  error "enable ENABLE_ATOM_TRACKING"
# endif
# if ENABLE_ECUT_CORRECTION
#  error "disable ENABLE_ECUT_CORRECTION"
# endif
# if ENABLE_NTOTATOMS_CORRECTION
#  error "disable ENABLE_NTOTATOMS_CORRECTION"
# endif
#endif

#if ENABLE_REF_ENERGY
# if !ENABLE_REALTIME_ENERGY
#  error "enable ENABLE_REALTIME_ENERGY"
# endif
#endif

#if ENABLE_REF_NEARBY
# if !ENABLE_REF_ATOMS
#  error "enable ENABLE_REF_ATOMS"
# endif
#endif

#endif /* CONFIG_H */
