static void add(void) {
  int i, j, k, m;

  // This loop nest performs element-wise addition: u[...] = u[...] + rhs[...]
  // Analyze for OpenMP suitability based on the objective:
  // 1. Data Dependencies:
  //    - Each iteration u[m][i][j][k] = u[m][i][j][k] + rhs[m][i][j][k]
  //      only depends on the values at the specific index (m, i, j, k).
  //    - There are NO loop-carried dependencies across any of the m, i, j, or k loops.
  //    - This means iterations are independent and highly suitable for parallel execution.
  //
  // 2. Load Imbalance:
  //    - Assuming grid_points dimensions and loop bounds are fixed and reasonable,
  //      the computation inside the inner loop (the addition) is constant time per element.
  //    - The workload for each iteration of any outer loop (m, i, j) is the sum
  //      of work in the inner loops, resulting in a relatively even distribution
  //      if parallelized across m, i, or j (or a combination).
  //    - No structural changes (like loop tiling or partitioning) are strictly needed
  //      for basic load balancing; standard OpenMP scheduling will work well.
  //
  // 3. Synchronization Overhead:
  //    - Each thread will be writing to a unique memory location u[m][i][j][k].
  //    - There are no shared writes or race conditions on 'u' or 'rhs'.
  //    - No explicit synchronization (locks, critical sections, atomics) is required
  //      for the addition operation itself.
  //    - The structure is already optimal regarding synchronization needs - there are none.
  //
  // Conclusion: The current loop structure is already highly parallelizable
  // without needing significant structural changes like dependency breaking
  // or complex synchronization handling. Parallelization can be applied directly
  // to one or more of the outer loops (m, i, j) using OpenMP directives like #pragma omp for.

  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = 1; k <= grid_points[2]-2; k++) {
	  // Independent operation per element
	  u[m][i][j][k] = u[m][i][j][k] + rhs[m][i][j][k];
	}
      }
    }
  }
}