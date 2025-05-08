static void z_solve_cell(void) {
  int i,j,k,ksize;
  ksize = grid_points[2]-1;
  #pragma omp parallel for private(j)
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      binvcrhs( lhs[i][j][0][BB],
		lhs[i][j][0][CC],
		rhs[i][j][0] );
    }
  }
  for (k = 1; k < ksize; k++) {
      for (i = 1; i < grid_points[0]-1; i++) {
	#pragma omp parallel for private(j)
	  for (j = 1; j < grid_points[1]-1; j++) {
	matvec_sub(lhs[i][j][k][AA],
		   rhs[i][j][k-1], rhs[i][j][k]);
	matmul_sub(lhs[i][j][k][AA],
		   lhs[i][j][k-1][CC],
		   lhs[i][j][k][BB]);
	binvcrhs( lhs[i][j][k][BB],
		  lhs[i][j][k][CC],
		  rhs[i][j][k] );
      }
    }
  }
  for (i = 1; i < grid_points[0]-1; i++) {
    #pragma omp parallel for private(j)
    for (j = 1; j < grid_points[1]-1; j++) {
      matvec_sub(lhs[i][j][ksize][AA],
		 rhs[i][j][ksize-1], rhs[i][j][ksize]);
      matmul_sub(lhs[i][j][ksize][AA],
		 lhs[i][j][ksize-1][CC],
		 lhs[i][j][ksize][BB]);
      binvrhs( lhs[i][j][ksize][BB],
	       rhs[i][j][ksize] );
    }
  }
}