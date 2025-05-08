static void y_backsubstitute(void) {
  int i, j, k, m, n;

  // The outermost j loop has a backward iteration with a loop-carried dependency:
  // rhs[i][j][k][m] depends on rhs[i][j+1][k][n]. This dependency prevents
  // parallelization of the j loop without a fundamental algorithmic change.
  // Therefore, the j loop remains sequential.
  for (j = grid_points[1]-2; j >= 0; j--) {

    // The loops i, k, and m iterate over independent target elements rhs[i][j][k][m]
    // for a fixed value of j. These loops represent the main computational workload
    // that can be distributed among threads.
    // OpenMP parallelization directives should target these loops, potentially
    // using collapse(3) or nested parallel loops depending on the desired granularity
    // and overhead considerations.
    // Example OpenMP Pragma Placement (Conceptual):
    // #pragma omp for collapse(3) private(i, k, m, n) // if inside parallel region
    // #pragma omp parallel for collapse(3) private(i, k, m, n) // if starting parallel region here
    for (i = 1; i < grid_points[0]-1; i++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	for (m = 0; m < BLOCK_SIZE; m++) {

          // The innermost n loop computes a sum that is subtracted from a single
          // element rhs[i][j][k][m]. Using a temporary variable (temp_sum) to
          // accumulate this sum makes the computation for a specific rhs[i][j][k][m]
          // self-contained. This prevents write conflicts if the i, k, m loops
          // are parallelized and avoids false sharing if temp_sum were part of rhs.
          // The read from rhs[i][j+1][k][n] accesses data from the previous (completed)
          // j iteration, which is safe in the sequential j loop context.
	  double temp_sum = 0.0;
	  for (n = 0; n < BLOCK_SIZE; n++) {
	    temp_sum += lhs[i][j][k][CC][m][n] * rhs[i][j+1][k][n];
	  }

          // The final update to rhs[i][j][k][m] uses the locally computed temp_sum.
          // This write is independent for distinct (i, k, m) tuples within the
          // same j iteration, aligning with the parallelization strategy for i, k, m.
	  rhs[i][j][k][m] = rhs[i][j][k][m] - temp_sum;
	}
      }
    }
  }
}