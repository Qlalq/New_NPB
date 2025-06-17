static void z_backsubstitute(void) {
  int i, j, k, m, n;
  
  #pragma omp parallel for private(i, j, k, m, n) collapse(2)
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = grid_points[2]-2; k >= 0; k--) {
        for (m = 0; m < BLOCK_SIZE; m++) {
          for (n = 0; n < BLOCK_SIZE; n++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] 
              - lhs[i][j][k][CC][m][n]*rhs[i][j][k+1][n];
          }
        }
      }
    }
  }
}