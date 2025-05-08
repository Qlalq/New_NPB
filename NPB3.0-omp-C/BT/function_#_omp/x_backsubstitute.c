static void x_backsubstitute(void) {
  int i, j, k, m, n;
  for (i = grid_points[0]-2; i >= 0; i--) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
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