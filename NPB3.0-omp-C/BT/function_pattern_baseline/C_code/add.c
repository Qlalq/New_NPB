static void add(void) {
  int i, j, k, m;
  
  #pragma omp parallel for private(j, k, m) collapse(3)
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
        for (m = 0; m < 5; m++) {
          u[i][j][k][m] = u[i][j][k][m] + rhs[i][j][k][m];
        }
      }
    }
  }
}