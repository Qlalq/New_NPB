static void y_solve_cell(void) {

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c     performs guaussian elimination on this cell.
c     
c     assumes that unpacking routines for non-first cells 
c     preload C' and rhs' from previous cell.
c     
c     assumed send happens outside this routine, but that
c     c'(JMAX) and rhs'(JMAX) will be sent to next cell
c-------------------------------------------------------------------*/

  int i, j, k, jsize;

  jsize = grid_points[1]-1;

#pragma omp for  
  for (i = 1; i < grid_points[0]-1; i++) {
    for (k = 1; k < grid_points[2]-1; k++) {

/*--------------------------------------------------------------------
c     multiply c(i,0,k) by b_inverse and copy back to c
c     multiply rhs(0) by b_inverse(0) and copy to rhs
c-------------------------------------------------------------------*/
      binvcrhs( lhs[i][0][k][BB],
		lhs[i][0][k][CC],
		rhs[i][0][k] );
    }
  }

/*--------------------------------------------------------------------
c     begin inner most do loop
c     do all the elements of the cell unless last 
c-------------------------------------------------------------------*/
  for (j = 1; j < jsize; j++) {
#pragma omp for
    for (i = 1; i < grid_points[0]-1; i++) {
      for (k = 1; k < grid_points[2]-1; k++) {

/*--------------------------------------------------------------------
c     subtract A*lhs_vector(j-1) from lhs_vector(j)
c     
c     rhs(j) = rhs(j) - A*rhs(j-1)
c-------------------------------------------------------------------*/
	matvec_sub(lhs[i][j][k][AA],
		   rhs[i][j-1][k], rhs[i][j][k]);

/*--------------------------------------------------------------------
c     B(j) = B(j) - C(j-1)*A(j)
c-------------------------------------------------------------------*/
	matmul_sub(lhs[i][j][k][AA],
		   lhs[i][j-1][k][CC],
		   lhs[i][j][k][BB]);

/*--------------------------------------------------------------------
c     multiply c(i,j,k) by b_inverse and copy back to c
c     multiply rhs(i,1,k) by b_inverse(i,1,k) and copy to rhs
c-------------------------------------------------------------------*/
	binvcrhs( lhs[i][j][k][BB],
		  lhs[i][j][k][CC],
		  rhs[i][j][k] );

      }
    }
  }

#pragma omp for  
  for (i = 1; i < grid_points[0]-1; i++) {
    for (k = 1; k < grid_points[2]-1; k++) {

/*--------------------------------------------------------------------
c     rhs(jsize) = rhs(jsize) - A*rhs(jsize-1)
c-------------------------------------------------------------------*/
      matvec_sub(lhs[i][jsize][k][AA],
		 rhs[i][jsize-1][k], rhs[i][jsize][k]);

/*--------------------------------------------------------------------
c     B(jsize) = B(jsize) - C(jsize-1)*A(jsize)
c     call matmul_sub(aa,i,jsize,k,c,
c     $              cc,i,jsize-1,k,c,BB,i,jsize,k)
c-------------------------------------------------------------------*/
      matmul_sub(lhs[i][jsize][k][AA],
		 lhs[i][jsize-1][k][CC],
		 lhs[i][jsize][k][BB]);

/*--------------------------------------------------------------------
c     multiply rhs(jsize) by b_inverse(jsize) and copy to rhs
c-------------------------------------------------------------------*/
      binvrhs( lhs[i][jsize][k][BB],
	       rhs[i][jsize][k] );

    }
  }
}