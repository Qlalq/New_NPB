static void z_solve(void) {

#pragma omp parallel
{

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c this function performs the solution of the approximate factorization
c step in the z-direction for all five matrix components
c simultaneously. The Thomas algorithm is employed to solve the
c systems for the z-lines. Boundary conditions are non-periodic
c-------------------------------------------------------------------*/

  int i, j, k, n, k1, k2, m;
  double fac1, fac2;

/*--------------------------------------------------------------------
c                          FORWARD ELIMINATION  
c-------------------------------------------------------------------*/

  lhsz();

  n = 0;

#pragma omp for
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 0; k <= grid_points[2]-3; k++) {
	k1 = k  + 1;
	k2 = k  + 2;
	fac1               = 1./lhs[n+2][i][j][k];
	lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k];
	lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k];
	for (m = 0; m < 3; m++) {
	  rhs[m][i][j][k] = fac1*rhs[m][i][j][k];
	}
	lhs[n+2][i][j][k1] = lhs[n+2][i][j][k1] -
	  lhs[n+1][i][j][k1]*lhs[n+3][i][j][k];
	lhs[n+3][i][j][k1] = lhs[n+3][i][j][k1] -
	  lhs[n+1][i][j][k1]*lhs[n+4][i][j][k];
	for (m = 0; m < 3; m++) {
	  rhs[m][i][j][k1] = rhs[m][i][j][k1] -
	    lhs[n+1][i][j][k1]*rhs[m][i][j][k];
	}
	lhs[n+1][i][j][k2] = lhs[n+1][i][j][k2] -
	  lhs[n+0][i][j][k2]*lhs[n+3][i][j][k];
	lhs[n+2][i][j][k2] = lhs[n+2][i][j][k2] -
	  lhs[n+0][i][j][k2]*lhs[n+4][i][j][k];
	for (m = 0; m < 3; m++) {
	  rhs[m][i][j][k2] = rhs[m][i][j][k2] -
	    lhs[n+0][i][j][k2]*rhs[m][i][j][k];
	}
      }
    }
  }

/*--------------------------------------------------------------------
c      The last two rows in this grid block are a bit different, 
c      since they do not have two more rows available for the
c      elimination of off-diagonal entries
c-------------------------------------------------------------------*/
  k  = grid_points[2]-2;
  k1 = grid_points[2]-1;
#pragma omp for
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      fac1               = 1./lhs[n+2][i][j][k];
      lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k];
      lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k];
      for (m = 0; m < 3; m++) {
	rhs[m][i][j][k] = fac1*rhs[m][i][j][k];
      }
      lhs[n+2][i][j][k1] = lhs[n+2][i][j][k1] -
	lhs[n+1][i][j][k1]*lhs[n+3][i][j][k];
      lhs[n+3][i][j][k1] = lhs[n+3][i][j][k1] -
	lhs[n+1][i][j][k1]*lhs[n+4][i][j][k];
      for (m = 0; m < 3; m++) {
	rhs[m][i][j][k1] = rhs[m][i][j][k1] -
	  lhs[n+1][i][j][k1]*rhs[m][i][j][k];
      }

/*--------------------------------------------------------------------
c               scale the last row immediately
c-------------------------------------------------------------------*/
      fac2               = 1./lhs[n+2][i][j][k1];
      for (m = 0; m < 3; m++) {
	rhs[m][i][j][k1] = fac2*rhs[m][i][j][k1];
      }
    }
  }

/*--------------------------------------------------------------------
c      do the u+c and the u-c factors               
c-------------------------------------------------------------------*/
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5;
#pragma omp for
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = 0; k <= grid_points[2]-3; k++) {
	k1 = k  + 1;
	k2 = k  + 2;
	  fac1               = 1./lhs[n+2][i][j][k];
	  lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k];
	  lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k];
	  rhs[m][i][j][k] = fac1*rhs[m][i][j][k];
	  lhs[n+2][i][j][k1] = lhs[n+2][i][j][k1] -
	    lhs[n+1][i][j][k1]*lhs[n+3][i][j][k];
	  lhs[n+3][i][j][k1] = lhs[n+3][i][j][k1] -
	    lhs[n+1][i][j][k1]*lhs[n+4][i][j][k];
	  rhs[m][i][j][k1] = rhs[m][i][j][k1] -
	    lhs[n+1][i][j][k1]*rhs[m][i][j][k];
	  lhs[n+1][i][j][k2] = lhs[n+1][i][j][k2] -
	    lhs[n+0][i][j][k2]*lhs[n+3][i][j][k];
	  lhs[n+2][i][j][k2] = lhs[n+2][i][j][k2] -
	    lhs[n+0][i][j][k2]*lhs[n+4][i][j][k];
	  rhs[m][i][j][k2] = rhs[m][i][j][k2] -
	    lhs[n+0][i][j][k2]*rhs[m][i][j][k];
	}
      }
    }

/*--------------------------------------------------------------------
c         And again the last two rows separately
c-------------------------------------------------------------------*/
    k  = grid_points[2]-2;
    k1 = grid_points[2]-1;
#pragma omp for
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	fac1               = 1./lhs[n+2][i][j][k];
	lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k];
	lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k];
	rhs[m][i][j][k]     = fac1*rhs[m][i][j][k];
	lhs[n+2][i][j][k1] = lhs[n+2][i][j][k1] -
	  lhs[n+1][i][j][k1]*lhs[n+3][i][j][k];
	lhs[n+3][i][j][k1] = lhs[n+3][i][j][k1] -
	  lhs[n+1][i][j][k1]*lhs[n+4][i][j][k];
	rhs[m][i][j][k1]   = rhs[m][i][j][k1] -
	  lhs[n+1][i][j][k1]*rhs[m][i][j][k];
/*--------------------------------------------------------------------
c               Scale the last row immediately (some of this is overkill
c               if this is the last cell)
c-------------------------------------------------------------------*/
	fac2               = 1./lhs[n+2][i][j][k1];
	rhs[m][i][j][k1]   = fac2*rhs[m][i][j][k1];

      }
    }
  }

/*--------------------------------------------------------------------
c                         BACKSUBSTITUTION 
c-------------------------------------------------------------------*/

  k  = grid_points[2]-2;
  k1 = grid_points[2]-1;
  n = 0;
  for (m = 0; m < 3; m++) {
#pragma omp for
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	rhs[m][i][j][k] = rhs[m][i][j][k] -
	  lhs[n+3][i][j][k]*rhs[m][i][j][k1];
      }
    }
  }

  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5;
#pragma omp for
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	rhs[m][i][j][k] = rhs[m][i][j][k] -
	  lhs[n+3][i][j][k]*rhs[m][i][j][k1];
      }
    }
  }

/*--------------------------------------------------------------------
c      Whether or not this is the last processor, we always have
c      to complete the back-substitution 
c-------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c      The first three factors
c-------------------------------------------------------------------*/
  n = 0;
  for (m = 0; m < 3; m++) {
#pragma omp for
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = grid_points[2]-3; k >= 0; k--) {
	  k1 = k  + 1;
	  k2 = k  + 2;
	  rhs[m][i][j][k] = rhs[m][i][j][k] - 
	    lhs[n+3][i][j][k]*rhs[m][i][j][k1] -
	    lhs[n+4][i][j][k]*rhs[m][i][j][k2];
	}
      }
    }
  }

/*--------------------------------------------------------------------
c      And the remaining two
c-------------------------------------------------------------------*/
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5;
#pragma omp for
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = grid_points[2]-3; k >= 0; k--) {
	  k1 = k  + 1;
	  k2 = k  + 2;
	  rhs[m][i][j][k] = rhs[m][i][j][k] - 
	    lhs[n+3][i][j][k]*rhs[m][i][j][k1] -
	    lhs[n+4][i][j][k]*rhs[m][i][j][k2];
	}
      }
    }
  }

}

  tzetar();
}