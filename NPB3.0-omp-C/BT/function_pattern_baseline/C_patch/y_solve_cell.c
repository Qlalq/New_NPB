static void y_solve_cell(void) {
  int i, j, k, jsize;

  jsize = grid_points[1]-1;
  #pragma omp parallel for private(i, k)
  for (i = 1; i < grid_points[0]-1; i++) {
    for (k = 1; k < grid_points[2]-1; k++) {
      binvcrhs( lhs[i][0][k][BB],
		lhs[i][0][k][CC],
		rhs[i][0][k] );
    }
  }
  for (j = 1; j < jsize; j++) {
    for (i = 1; i < grid_points[0]-1; i++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	matvec_sub(lhs[i][j][k][AA],
		   rhs[i][j-1][k], rhs[i][j][k]);
	matmul_sub(lhs[i][j][k][AA],
		   lhs[i][j-1][k][CC],
		   lhs[i][j][k][BB]);
	binvcrhs( lhs[i][j][k][BB],
		  lhs[i][j][k][CC],
		  rhs[i][j][k] );
      }
    }
  }
  for (i = 1; i < grid_points[0]-1; i++) {
    #pragma omp parallel for private(i, k)
    for (k = 1; k < grid_points[2]-1; k++) {
      matvec_sub(lhs[i][jsize][k][AA],
		 rhs[i][jsize-1][k], rhs[i][jsize][k]);
      matmul_sub(lhs[i][jsize][k][AA],
		 lhs[i][jsize-1][k][CC],
		 lhs[i][jsize][k][BB]);
      binvrhs( lhs[i][jsize][k][BB],
	       rhs[i][jsize][k] );
    }
  }
}