static void x_solve(void) {

#pragma omp parallel
{

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c this function performs the solution of the approximate factorization
c step in the x-direction for all five matrix components
c simultaneously. The Thomas algorithm is employed to solve the
c systems for the x-lines. Boundary conditions are non-periodic
--------------------------------------------------------------------*/

  int i, j, k, n, i1, i2, m;
  double fac1, fac2;

/*--------------------------------------------------------------------
c                          FORWARD ELIMINATION  
--------------------------------------------------------------------*/
  lhsx();

/*--------------------------------------------------------------------
c      perform the Thomas algorithm; first, FORWARD ELIMINATION     
--------------------------------------------------------------------*/
  n = 0;
  for (i = 0; i <= grid_points[0]-3; i++) {
    i1 = i  + 1;
    i2 = i  + 2;
#pragma omp for
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	fac1               = 1./lhs[n+2][i][j][k];
	lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k];
	lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k];
	for (m = 0; m < 3; m++) {
	  rhs[m][i][j][k] = fac1*rhs[m][i][j][k];
	}
	lhs[n+2][i1][j][k] = lhs[n+2][i1][j][k] -
	  lhs[n+1][i1][j][k]*lhs[n+3][i][j][k];
	lhs[n+3][i1][j][k] = lhs[n+3][i1][j][k] -
	  lhs[n+1][i1][j][k]*lhs[n+4][i][j][k];
	for (m = 0; m < 3; m++) {
	  rhs[m][i1][j][k] = rhs[m][i1][j][k] -
	    lhs[n+1][i1][j][k]*rhs[m][i][j][k];
	}
	lhs[n+1][i2][j][k] = lhs[n+1][i2][j][k] -
	  lhs[n+0][i2][j][k]*lhs[n+3][i][j][k];
	lhs[n+2][i2][j][k] = lhs[n+2][i2][j][k] -
	  lhs[n+0][i2][j][k]*lhs[n+4][i][j][k];
	for (m = 0; m < 3; m++) {
	  rhs[m][i2][j][k] = rhs[m][i2][j][k] -
	    lhs[n+0][i2][j][k]*rhs[m][i][j][k];
	}
      }
    }
  }

/*--------------------------------------------------------------------
c      The last two rows in this grid block are a bit different, 
c      since they do not have two more rows available for the
c      elimination of off-diagonal entries
--------------------------------------------------------------------*/

  i  = grid_points[0]-2;
  i1 = grid_points[0]-1;
#pragma omp for
  for (j = 1; j <= grid_points[1]-2; j++) {
    for (k = 1; k <= grid_points[2]-2; k++) {
      fac1               = 1.0/lhs[n+2][i][j][k];
      lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k];
      lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k];
      for (m = 0; m < 3; m++) {
	rhs[m][i][j][k] = fac1*rhs[m][i][j][k];
      }
      lhs[n+2][i1][j][k] = lhs[n+2][i1][j][k] -
	lhs[n+1][i1][j][k]*lhs[n+3][i][j][k];
      lhs[n+3][i1][j][k] = lhs[n+3][i1][j][k] -
	lhs[n+1][i1][j][k]*lhs[n+4][i][j][k];
      for (m = 0; m < 3; m++) {
	rhs[m][i1][j][k] = rhs[m][i1][j][k] -
	  lhs[n+1][i1][j][k]*rhs[m][i][j][k];
      }

/*--------------------------------------------------------------------
c            scale the last row immediately 
--------------------------------------------------------------------*/
      fac2               = 1./lhs[n+2][i1][j][k];
      for (m = 0; m < 3; m++) {
	rhs[m][i1][j][k] = fac2*rhs[m][i1][j][k];
      }
    }
  }

/*--------------------------------------------------------------------
c      do the u+c and the u-c factors                 
--------------------------------------------------------------------*/

  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5;
    for (i = 0; i <= grid_points[0]-3; i++) {
      i1 = i  + 1;
      i2 = i  + 2;
#pragma omp for
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = 1; k <= grid_points[2]-2; k++) {
	  fac1               = 1./lhs[n+2][i][j][k];
	  lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k];
	  lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k];
	  rhs[m][i][j][k] = fac1*rhs[m][i][j][k];
	  lhs[n+2][i1][j][k] = lhs[n+2][i1][j][k] -
	    lhs[n+1][i1][j][k]*lhs[n+3][i][j][k];
	  lhs[n+3][i1][j][k] = lhs[n+3][i1][j][k] -
	    lhs[n+1][i1][j][k]*lhs[n+4][i][j][k];
	  rhs[m][i1][j][k] = rhs[m][i1][j][k] -
	    lhs[n+1][i1][j][k]*rhs[m][i][j][k];
	  lhs[n+1][i2][j][k] = lhs[n+1][i2][j][k] -
	    lhs[n+0][i2][j][k]*lhs[n+3][i][j][k];
	  lhs[n+2][i2][j][k] = lhs[n+2][i2][j][k] -
	    lhs[n+0][i2][j][k]*lhs[n+4][i][j][k];
	  rhs[m][i2][j][k] = rhs[m][i2][j][k] -
	    lhs[n+0][i2][j][k]*rhs[m][i][j][k];
	}
      }
    }

/*--------------------------------------------------------------------
c         And again the last two rows separately
--------------------------------------------------------------------*/
    i  = grid_points[0]-2;
    i1 = grid_points[0]-1;
    
#pragma omp for
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	fac1               = 1./lhs[n+2][i][j][k];
	lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k];
	lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k];
	rhs[m][i][j][k]     = fac1*rhs[m][i][j][k];
	lhs[n+2][i1][j][k] = lhs[n+2][i1][j][k] -
	  lhs[n+1][i1][j][k]*lhs[n+3][i][j][k];
	lhs[n+3][i1][j][k] = lhs[n+3][i1][j][k] -
	  lhs[n+1][i1][j][k]*lhs[n+4][i][j][k];
	rhs[m][i1][j][k]   = rhs[m][i1][j][k] -
	  lhs[n+1][i1][j][k]*rhs[m][i][j][k];
/*--------------------------------------------------------------------
c               Scale the last row immediately
--------------------------------------------------------------------*/
	fac2               = 1./lhs[n+2][i1][j][k];
	rhs[m][i1][j][k]   = fac2*rhs[m][i1][j][k];

      }
    }
  }

/*--------------------------------------------------------------------
c                         BACKSUBSTITUTION 
--------------------------------------------------------------------*/

  i  = grid_points[0]-2;
  i1 = grid_points[0]-1;
  n = 0;
  for (m = 0; m < 3; m++) {
#pragma omp for
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	rhs[m][i][j][k] = rhs[m][i][j][k] -
	  lhs[n+3][i][j][k]*rhs[m][i1][j][k];
      }
    }
  }

  for (m = 3; m < 5; m++) {
#pragma omp for
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	n = (m-3+1)*5;
	rhs[m][i][j][k] = rhs[m][i][j][k] -
	  lhs[n+3][i][j][k]*rhs[m][i1][j][k];
      }
    }
  }

/*--------------------------------------------------------------------
c      The first three factors
--------------------------------------------------------------------*/
  n = 0;
  for (i = grid_points[0]-3; i >= 0; i--) {
    i1 = i  + 1;
    i2 = i  + 2;
#pragma omp for
    for (m = 0; m < 3; m++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = 1; k <= grid_points[2]-2; k++) {
	  rhs[m][i][j][k] = rhs[m][i][j][k] - 
	    lhs[n+3][i][j][k]*rhs[m][i1][j][k] -
	    lhs[n+4][i][j][k]*rhs[m][i2][j][k];
	}
      }
    }
  }

/*--------------------------------------------------------------------
c      And the remaining two
--------------------------------------------------------------------*/
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5;
    for (i = grid_points[0]-3; i >= 0; i--) {
      i1 = i  + 1;
      i2 = i  + 2;
#pragma omp for
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = 1; k <= grid_points[2]-2; k++) {
	  rhs[m][i][j][k] = rhs[m][i][j][k] - 
	    lhs[n+3][i][j][k]*rhs[m][i1][j][k] -
	    lhs[n+4][i][j][k]*rhs[m][i2][j][k];
	}
      }
    }
  }

}

/*--------------------------------------------------------------------
c      Do the block-diagonal inversion          
--------------------------------------------------------------------*/
  ninvr();
}