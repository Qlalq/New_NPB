static void evolve(dcomplex u0[NZ][NY][NX], dcomplex u1[NZ][NY][NX],
                   int t, int indexmap[NZ][NY][NX], int d[3]) {
    int i, j, k;
    // Using explicit local variables for dimension sizes improves clarity
    // and makes the loop bounds explicit for parallelization directives.
    int nx = d[0]; // Size of the innermost loop dimension (e.g., X)
    int ny = d[1]; // Size of the middle loop dimension (e.g., Y)
    int nz = d[2]; // Size of the outermost loop dimension (e.g., Z)

    // This triple nested loop iterates over a 3D grid.
    // Inside the innermost loop, the operation is:
    // u1[k][j][i] = u0[k][j][i] * ex[t * indexmap[k][j][i]] (conceptually, via crmul)
    //
    // 1. Data Dependencies:
    //    - The write operation is to u1[k][j][i]. For each unique combination of (k, j, i),
    //      a unique memory location in u1 is written to.
    //    - The read operations are from u0[k][j][i] and ex[...]. These are read-only for
    //      the purpose of computing u1[k][j][i].
    //    - There are no loop-carried dependencies. The computation for u1[k][j][i] does not
    //      depend on the value of u1 from a previous iteration (e.g., u1[k][j][i-1] or u1[k-1][j][i]).
    //    - The index for ex depends on t and indexmap[k][j][i]. While indexmap is read,
    //      its value for a given (k, j, i) is constant for that iteration and does not create
    //      a loop-carried dependency.
    //    - Therefore, each iteration (k, j, i) is independent and can be executed in any order.
    //      No structural change is needed to break dependencies.
    //
    // 2. Load Imbalance:
    //    - Assuming the cost of the crmul operation is roughly constant for all valid inputs,
    //      the workload for each iteration (k, j, i) is approximately uniform.
    //    - Standard OpenMP loop parallelization directives (like #pragma omp parallel for)
    //      with default scheduling or static/dynamic scheduling will distribute the iterations
    //      across threads, providing reasonable load balance for a uniform workload.
    //    - No structural change like loop partitioning or adaptive scheduling is needed at
    //      this level of refinement, as the uniform workload is handled well by OpenMP.
    //
    // 3. Synchronization Overhead:
    //    - Since each iteration writes to a unique memory location u1[k][j][i], there are
    //      no concurrent writes to shared data by different iterations.
    //    - No explicit synchronization (like locks, atomics, critical sections) is required
    //      within the loop body.
    //    - The structure inherently minimizes synchronization needs, requiring only potential
    //      implicit synchronization at the end of the parallel region (e.g., via a barrier).
    //
    // Conclusion: The original loop structure is already well-suited for OpenMP parallelization
    // because the iterations are independent and the workload is uniform. The primary refinement
    // applied here is making the loop bounds clearer using local variables, which aids in
    // writing clear OpenMP directives (e.g., for collapse clauses). No fundamental structural
    // changes (like introducing temporary arrays to break dependencies) are necessary as
    // no such dependencies exist.

    for (k = 0; k < nz; k++) {
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                crmul(u1[k][j][i], u0[k][j][i], ex[t*indexmap[k][j][i]]);
            }
        }
    }
}