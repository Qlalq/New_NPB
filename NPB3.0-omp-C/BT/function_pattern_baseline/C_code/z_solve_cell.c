static void z_solve_cell(void) {
  int i, j, k, ksize;
  ksize = grid_points[2] - 1;

  /* solve at k = 0 */
  #pragma omp for collapse(2) schedule(static)
  for (i = 1; i < grid_points[0] - 1; i++) {
    for (j = 1; j < grid_points[1] - 1; j++) {
      binvcrhs(lhs[i][j][0][BB],
               lhs[i][j][0][CC],
               rhs[i][j][0]);
    }
  }

  /* forward elimination for k = 1 … ksize-1 */
  for (k = 1; k < ksize; k++) {
    #pragma omp for collapse(2) schedule(static)
    for (i = 1; i < grid_points[0] - 1; i++) {
      for (j = 1; j < grid_points[1] - 1; j++) {
        matvec_sub(lhs[i][j][k][AA],
                   rhs[i][j][k - 1],
                   rhs[i][j][k]);
        matmul_sub(lhs[i][j][k][AA],
                   lhs[i][j][k - 1][CC],
                   lhs[i][j][k][BB]);
        binvcrhs(lhs[i][j][k][BB],
                 lhs[i][j][k][CC],
                 rhs[i][j][k]);
      }
    }
    /* implicit barrier here ensures all threads finish k before next */
  }

  /* solve at k = ksize */
  #pragma omp for collapse(2) schedule(static)
  for (i = 1; i < grid_points[0] - 1; i++) {
    for (j = 1; j < grid_points[1] - 1; j++) {
      matvec_sub(lhs[i][j][ksize][AA],
                 rhs[i][j][ksize - 1],
                 rhs[i][j][ksize]);
      matmul_sub(lhs[i][j][ksize][AA],
                 lhs[i][j][ksize - 1][CC],
                 lhs[i][j][ksize][BB]);
      binvrhs(lhs[i][j][ksize][BB],
              rhs[i][j][ksize]);
    }
  }
}