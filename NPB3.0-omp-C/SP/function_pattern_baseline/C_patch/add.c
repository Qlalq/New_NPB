static void add(void) {
  int i, j, k, m;
#pragma omp parallel for collapse(4)
  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = 1; k <= grid_points[2]-2; k++) {
	  u[m][i][j][k] = u[m][i][j][k] + rhs[m][i][j][k];
	}
      }
    }
  }
}