static void z_backsubstitute(void) {
  int i, j, k, m, n;

  // The loops over i and j are independent and good candidates for
  // parallelization using OpenMP directives like #pragma omp parallel for
  // or #pragma omp tasks.
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {

      // The loop over k proceeds backward and exhibits a loop-carried
      // dependency (computation for k depends on results from k+1).
      // This loop is generally sequential or requires specialized parallel
      // algorithms for backward substitution (e.g., wavefront, recursive doubling).
      for (k = grid_points[2]-2; k >= 0; k--) {

        // For a fixed (i, j, k), the loop over m iterates over
        // independent elements of rhs[i][j][k][:]. This loop is
        // parallelizable (e.g., #pragma omp parallel for).
	for (m = 0; m < BLOCK_SIZE; m++) {

          // The innermost loop over n computes a sum of terms that are
          // subtracted from rhs[i][j][k][m]. This is a reduction-like
          // operation.
          // Using a temporary variable makes the reduction structure explicit.
          double sum_of_terms = 0.0; // Use appropriate floating-point type
          // This loop can be parallelized as a reduction (#pragma omp parallel for reduction(+:sum_of_terms))
	  for (n = 0; n < BLOCK_SIZE; n++) {
            // Assuming lhs and rhs hold floating-point values
	    sum_of_terms += lhs[i][j][k][CC][m][n]*rhs[i][j][k+1][n];
	  }

          // Apply the accumulated sum to the result element.
	  rhs[i][j][k][m] -= sum_of_terms;
	}
      }
    }
  }
}