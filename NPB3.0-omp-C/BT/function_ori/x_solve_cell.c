static void x_solve_cell(void) {

/*--------------------------------------------------------------------
c     performs guaussian elimination on this cell.
c     
c     assumes that unpacking routines for non-first cells 
c     preload C' and rhs' from previous cell.
c     
c     assumed send happens outside this routine, but that
c     c'(IMAX) and rhs'(IMAX) will be sent to next cell
c-------------------------------------------------------------------*/

  int i,j,k,isize;

  isize = grid_points[0]-1;

/*--------------------------------------------------------------------
c     outer most do loops - sweeping in i direction
c-------------------------------------------------------------------*/
#pragma omp for
  for (j = 1; j < grid_points[1]-1; j++) {
    for (k = 1; k < grid_points[2]-1; k++) {

/*--------------------------------------------------------------------
c     multiply c(0,j,k) by b_inverse and copy back to c
c     multiply rhs(0) by b_inverse(0) and copy to rhs
c-------------------------------------------------------------------*/
      binvcrhs( lhs[0][j][k][BB],
		lhs[0][j][k][CC],
		rhs[0][j][k] );
    }
  }

/*--------------------------------------------------------------------
c     begin inner most do loop
c     do all the elements of the cell unless last 
c-------------------------------------------------------------------*/
  for (i = 1; i < isize; i++) {
#pragma omp for
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {

/*--------------------------------------------------------------------
c     rhs(i) = rhs(i) - A*rhs(i-1)
c-------------------------------------------------------------------*/
	matvec_sub(lhs[i][j][k][AA],
		   rhs[i-1][j][k], rhs[i][j][k]);

/*--------------------------------------------------------------------
c     B(i) = B(i) - C(i-1)*A(i)
c-------------------------------------------------------------------*/
	matmul_sub(lhs[i][j][k][AA],
		   lhs[i-1][j][k][CC],
		   lhs[i][j][k][BB]);


/*--------------------------------------------------------------------
c     multiply c(i,j,k) by b_inverse and copy back to c
c     multiply rhs(1,j,k) by b_inverse(1,j,k) and copy to rhs
c-------------------------------------------------------------------*/
	binvcrhs( lhs[i][j][k][BB],
		  lhs[i][j][k][CC],
		  rhs[i][j][k] );

      }
    }
  }

#pragma omp for
  for (j = 1; j < grid_points[1]-1; j++) {
    for (k = 1; k < grid_points[2]-1; k++) {

/*--------------------------------------------------------------------
c     rhs(isize) = rhs(isize) - A*rhs(isize-1)
c-------------------------------------------------------------------*/
      matvec_sub(lhs[isize][j][k][AA],
		 rhs[isize-1][j][k], rhs[isize][j][k]);

/*--------------------------------------------------------------------
c     B(isize) = B(isize) - C(isize-1)*A(isize)
c-------------------------------------------------------------------*/
      matmul_sub(lhs[isize][j][k][AA],
		 lhs[isize-1][j][k][CC],
		 lhs[isize][j][k][BB]);

/*--------------------------------------------------------------------
c     multiply rhs() by b_inverse() and copy to rhs
c-------------------------------------------------------------------*/
      binvrhs( lhs[i][j][k][BB],
	       rhs[i][j][k] );

    }
  }
}