static void lhsinit(void) {
  int i, j, k, n;
  for (n = 0; n < 15; n++) {
    for (i = 0; i < grid_points[0]; i++) {
      for (j = 0; j < grid_points[1]; j++) {
	for (k = 0; k < grid_points[2]; k++) {
	  lhs[n][i][j][k] = 0.0;
	}
      }
    }
  }
  for (n = 0; n < 3; n++) {
    for (i = 0; i < grid_points[0]; i++) {
      for (j = 0; j < grid_points[1]; j++) {
	for (k = 0; k < grid_points[2]; k++) {
	  lhs[5*n+2][i][j][k] = 1.0;
	}
      }
    }
  }
}