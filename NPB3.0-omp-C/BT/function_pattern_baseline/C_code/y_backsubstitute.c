static void y_backsubstitute(void) {
  int i, j, k, m, n;
  
  #pragma omp parallel for private(i, j, k, m, n) collapse(3)
  for (j = grid_points[1]-2; j >= 0; j--) {
    for (i = 1; i < grid_points[0]-1; i++) {
      for (k = 1; k < grid_points[2]-1; k++) {
        for (m = 0; m < BLOCK_SIZE; m++) {
          for (n = 0; n < BLOCK_SIZE; n++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] 
              - lhs[i][j][k][CC][m][n] * rhs[i][j+1][k][n];
          }
        }
      }
    }
  }
}