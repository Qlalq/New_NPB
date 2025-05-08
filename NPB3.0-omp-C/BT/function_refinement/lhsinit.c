static void lhsinit(void) {
  int i, j, k, m, n;

  // Initialize lhs elements to 0.0.
  // Each iteration of the outer loops (i, j, k) operates on a distinct block
  // of the lhs array (lhs[i][j][k][...]).
  // This structure is highly parallelizable over i, j, or k.
  // No loop-carried dependencies within these loops.
  for (i = 0; i < grid_points[0]; i++) {
    for (j = 0; j < grid_points[1]; j++) {
      for (k = 0; k < grid_points[2]; k++) {
        for (m = 0; m < 5; m++) {
          for (n = 0; n < 5; n++) {
            lhs[i][j][k][0][m][n] = 0.0;
            lhs[i][j][k][1][m][n] = 0.0;
            lhs[i][j][k][2][m][n] = 0.0;
          }
        }
      }
    }
  }

  // Set diagonal elements of the second component (index 1) to 1.0.
  // Similar to the first block, each iteration of the outer loops (i, j, k)
  // operates on a distinct block.
  // This block is highly parallelizable over i, j, or k.
  // No loop-carried dependencies within these loops.
  // This block must execute after the first block completes the zero-initialization.
  for (i = 0; i < grid_points[0]; i++) {
    for (j = 0; j < grid_points[1]; j++) {
      for (k = 0; k < grid_points[2]; k++) {
        for (m = 0; m < 5; m++) {
          lhs[i][j][k][1][m][m] = 1.0;
        }
      }
    }
  }
}