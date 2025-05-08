static void x_backsubstitute(void) {
  int i, j, k, m, n;
  // The outermost loop over 'i' has a loop-carried dependency:
  // rhs[i] depends on rhs[i+1]. This loop must be executed sequentially.
  for (i = grid_points[0]-2; i >= 0; i--) {
    // The loops over 'j' and 'k' iterate over spatial dimensions.
    // For a fixed 'i', the computations for different (j, k) pairs
    // update independent elements of rhs[i].
    // These loops are parallelizable. They can be parallelized separately
    // or collapsed together for better load balancing depending on grid sizes.
    // OpenMP directive would be placed here, potentially with collapse(2).
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
        // The inner loops over 'm' and 'n' perform the core computation
        // for a specific rhs[i][j][k][m] element. The update is
        // a reduction-like sum over 'n' for each 'm'.
        // These loops are small and tightly coupled to the computation
        // for a single (i, j, k) point, best kept sequential within
        // the parallel task for (j, k).
	for (m = 0; m < BLOCK_SIZE; m++) {
	  for (n = 0; n < BLOCK_SIZE; n++) {
	    rhs[i][j][k][m] = rhs[i][j][k][m]
	      - lhs[i][j][k][CC][m][n]*rhs[i+1][j][k][n];
	  }
	}
      }
    }
  }
}