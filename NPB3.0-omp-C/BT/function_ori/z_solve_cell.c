static void z_solve_cell(void) {

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c     performs guaussian elimination on this cell.
c     
c     assumes that unpacking routines for non-first cells 
c     preload C' and rhs' from previous cell.
c     
c     assumed send happens outside this routine, but that
c     c'(KMAX) and rhs'(KMAX) will be sent to next cell.
c-------------------------------------------------------------------*/

  int i,j,k,ksize;

  ksize = grid_points[2]-1;

/*--------------------------------------------------------------------
c     outer most do loops - sweeping in i direction
c-------------------------------------------------------------------*/
#pragma omp for  
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {

/*--------------------------------------------------------------------
c     multiply c(i,j,0) by b_inverse and copy back to c
c     multiply rhs(0) by b_inverse(0) and copy to rhs
c-------------------------------------------------------------------*/
      binvcrhs( lhs[i][j][0][BB],
		lhs[i][j][0][CC],
		rhs[i][j][0] );

    }
  }

/*--------------------------------------------------------------------
c     begin inner most do loop
c     do all the elements of the cell unless last 
c-------------------------------------------------------------------*/
  for (k = 1; k < ksize; k++) {
#pragma omp for
      for (i = 1; i < grid_points[0]-1; i++) {
	  for (j = 1; j < grid_points[1]-1; j++) {

/*--------------------------------------------------------------------
c     subtract A*lhs_vector(k-1) from lhs_vector(k)
c     
c     rhs(k) = rhs(k) - A*rhs(k-1)
c-------------------------------------------------------------------*/
	matvec_sub(lhs[i][j][k][AA],
		   rhs[i][j][k-1], rhs[i][j][k]);

/*--------------------------------------------------------------------
c     B(k) = B(k) - C(k-1)*A(k)
c     call matmul_sub(aa,i,j,k,c,cc,i,j,k-1,c,BB,i,j,k)
c-------------------------------------------------------------------*/
	matmul_sub(lhs[i][j][k][AA],
		   lhs[i][j][k-1][CC],
		   lhs[i][j][k][BB]);

/*--------------------------------------------------------------------
c     multiply c(i,j,k) by b_inverse and copy back to c
c     multiply rhs(i,j,1) by b_inverse(i,j,1) and copy to rhs
c-------------------------------------------------------------------*/
	binvcrhs( lhs[i][j][k][BB],
		  lhs[i][j][k][CC],
		  rhs[i][j][k] );

      }
    }
  }

/*--------------------------------------------------------------------
c     Now finish up special cases for last cell
c-------------------------------------------------------------------*/
#pragma omp for  
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {

/*--------------------------------------------------------------------
c     rhs(ksize) = rhs(ksize) - A*rhs(ksize-1)
c-------------------------------------------------------------------*/
      matvec_sub(lhs[i][j][ksize][AA],
		 rhs[i][j][ksize-1], rhs[i][j][ksize]);

/*--------------------------------------------------------------------
c     B(ksize) = B(ksize) - C(ksize-1)*A(ksize)
c     call matmul_sub(aa,i,j,ksize,c,
c     $              cc,i,j,ksize-1,c,BB,i,j,ksize)
c-------------------------------------------------------------------*/
      matmul_sub(lhs[i][j][ksize][AA],
		 lhs[i][j][ksize-1][CC],
		 lhs[i][j][ksize][BB]);

/*--------------------------------------------------------------------
c     multiply rhs(ksize) by b_inverse(ksize) and copy to rhs
c-------------------------------------------------------------------*/
      binvrhs( lhs[i][j][ksize][BB],
	       rhs[i][j][ksize] );

    }
  }
}