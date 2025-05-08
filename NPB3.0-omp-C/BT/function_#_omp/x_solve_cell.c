static void x_solve_cell(void) {
  int i,j,k,isize;
  isize = grid_points[0]-1;
  for (j = 1; j < grid_points[1]-1; j++) {
    for (k = 1; k < grid_points[2]-1; k++) {
      binvcrhs( lhs[0][j][k][BB],
		lhs[0][j][k][CC],
		rhs[0][j][k] );
    }
  }
  for (i = 1; i < isize; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	matvec_sub(lhs[i][j][k][AA],
		   rhs[i-1][j][k], rhs[i][j][k]);
	matmul_sub(lhs[i][j][k][AA],
		   lhs[i-1][j][k][CC],
		   lhs[i][j][k][BB]);
	binvcrhs( lhs[i][j][k][BB],
		  lhs[i][j][k][CC],
		  rhs[i][j][k] );
      }
    }
  }
  for (j = 1; j < grid_points[1]-1; j++) {
    for (k = 1; k < grid_points[2]-1; k++) {
      matvec_sub(lhs[isize][j][k][AA],
		 rhs[isize-1][j][k], rhs[isize][j][k]);
      matmul_sub(lhs[isize][j][k][AA],
		 lhs[isize-1][j][k][CC],
		 lhs[isize][j][k][BB]);
      binvrhs( lhs[i][j][k][BB],
	       rhs[i][j][k] );
    }
  }
}