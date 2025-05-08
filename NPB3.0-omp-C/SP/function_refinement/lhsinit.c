static void lhsinit(void) {
  int i, j, k, n;

  // Phase 1: Initialize all relevant elements to 0.0
  // This loop iterates over 15 slices (n=0..14) of the LHS array.
  // Each slice is a 3D grid.
  // The operations within the inner loops (i, j, k) for a fixed n are independent.
  // The operations for different values of n are also independent.
  // This block is suitable for parallelizing the outermost 'n' loop
  // or the spatial loops (i, j, k), or a combination (e.g., collapse).
  for (n = 0; n < 15; n++) {
    for (i = 0; i < grid_points[0]; i++) {
      for (j = 0; j < grid_points[1]; j++) {
        for (k = 0; k < grid_points[2]; k++) {
          lhs[n][i][j][k] = 0.0;
        }
      }
    }
  }

  // Phase 2: Overwrite specific elements with 1.0
  // This loop iterates over 3 slices (n=0..2, corresponding to indices 2, 7, 12).
  // Each slice is a 3D grid.
  // This phase *must* happen after Phase 1 for the indices 2, 7, 12.
  // The operations within the inner loops (i, j, k) for a fixed n are independent.
  // The operations for different values of n (corresponding to indices 2, 7, 12)
  // are independent of each other.
  // This block is suitable for parallelizing the outermost 'n' loop
  // or the spatial loops (i, j, k), or a combination.
  for (n = 0; n < 3; n++) {
    // Calculate the target index in the first dimension (n)
    int target_n = 5 * n + 2; // This will be 2, 7, 12

    for (i = 0; i < grid_points[0]; i++) {
      for (j = 0; j < grid_points[1]; j++) {
        for (k = 0; k < grid_points[2]; k++) {
          // Write to the specific slice targeted by target_n
          // This write is independent for different (i, j, k) and different n (2, 7, 12).
          lhs[target_n][i][j][k] = 1.0;
        }
      }
    }
  }
}