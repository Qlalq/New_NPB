#include <stdio.h>

// Forward declarations for the functions called within adi
// (Assuming these functions are defined elsewhere and operate on shared data,
// which is common for ADI implementations like those in NPB)
extern void compute_rhs(void);
extern void txinvr(void);
extern void x_solve(void);
extern void y_solve(void);
extern void z_solve(void);
extern void add(void);

/**
 * @brief Performs one Alternating Direction Implicit (ADI) time step.
 *
 * This function orchestrates the sequence of kernel operations required
 * for a single ADI step. The steps are inherently sequential due to
 * data dependencies: each solve step (x, y, z) depends on the results
 * of the previous solve step applied to the computational state.
 *
 * Parallelism is typically applied *within* each of the called functions
 * (compute_rhs, txinvr, x_solve, y_solve, z_solve, add) by parallelizing
 * their internal loops or operations across independent data partitions
 * (e.g., parallelizing solves along one dimension by distributing loops
 * over the other two dimensions).
 *
 * This function is structured as a clear sequence of calls, making the
 * overall algorithm flow evident and guiding where internal parallelism
 * (e.g., using OpenMP directives within compute_rhs, x_solve loops, etc.)
 * should be introduced.
 */
static void adi(void) {
    // Step 1: Compute the right-hand side (RHS) of the implicit equations.
    // This step can often be parallelized by distributing the computation
    // across the grid points or computational domain.
    compute_rhs();

    // Step 2: Compute the inverse of the tri-diagonal matrices or factors
    // needed for the subsequent solve steps. Parallelism opportunities
    // depend on the specific factorization method used.
    txinvr();

    // Step 3: Perform the implicit solve along the X-direction.
    // This step involves solving many independent 1D tri-diagonal systems
    // along the X-lines. Parallelism is achieved by distributing the loops
    // over the Y and Z dimensions, allowing multiple X-solves to run concurrently.
    // This step depends on the output of compute_rhs and txinvr.
    x_solve();

    // Step 4: Perform the implicit solve along the Y-direction.
    // Solves along Y-lines based on the updated state from the X-solve.
    // Parallelism is achieved by distributing the loops over the X and Z dimensions.
    // This step depends on the output of x_solve.
    y_solve();

    // Step 5: Perform the implicit solve along the Z-direction.
    // Solves along Z-lines based on the updated state from the Y-solve.
    // Parallelism is achieved by distributing the loops over the X and Y dimensions.
    // This step depends on the output of y_solve.
    z_solve();

    // Step 6: Add the corrections or updates computed by the solve steps
    // to the main computational state. This step often involves simple
    // element-wise operations and can be parallelized by distributing
    // the loops over the grid points.
    // This step depends on the output of z_solve.
    add();
}