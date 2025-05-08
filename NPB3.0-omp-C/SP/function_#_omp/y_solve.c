static void y_solve(void) {
{
  int i, j, k, n, j1, j2, m;
  double fac1, fac2;
  lhsy();
  n = 0;
  for (j = 0; j <= grid_points[1]-3; j++) {
    j1 = j  + 1;
    j2 = j  + 2;
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	fac1               = 1./lhs[n+2][i][j][k];
	lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k];
	lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k];
        for (m = 0; m < 3; m++) {
	  rhs[m][i][j][k] = fac1*rhs[m][i][j][k];
	}
	lhs[n+2][i][j1][k] = lhs[n+2][i][j1][k] -
	  lhs[n+1][i][j1][k]*lhs[n+3][i][j][k];
	lhs[n+3][i][j1][k] = lhs[n+3][i][j1][k] -
	  lhs[n+1][i][j1][k]*lhs[n+4][i][j][k];
	for (m = 0; m < 3; m++) {
	  rhs[m][i][j1][k] = rhs[m][i][j1][k] -
	    lhs[n+1][i][j1][k]*rhs[m][i][j][k];
	}
	lhs[n+1][i][j2][k] = lhs[n+1][i][j2][k] -
	  lhs[n+0][i][j2][k]*lhs[n+3][i][j][k];
	lhs[n+2][i][j2][k] = lhs[n+2][i][j2][k] -
	  lhs[n+0][i][j2][k]*lhs[n+4][i][j][k];
	for (m = 0; m < 3; m++) {
	  rhs[m][i][j2][k] = rhs[m][i][j2][k] -
	    lhs[n+0][i][j2][k]*rhs[m][i][j][k];
	}
      }
    }
  }
  j  = grid_points[1]-2;
  j1 = grid_points[1]-1;
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (k = 1; k <= grid_points[2]-2; k++) {
      fac1               = 1./lhs[n+2][i][j][k];
      lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k];
      lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k];
      for (m = 0; m < 3; m++) {
	rhs[m][i][j][k] = fac1*rhs[m][i][j][k];
      }
      lhs[n+2][i][j1][k] = lhs[n+2][i][j1][k] -
	lhs[n+1][i][j1][k]*lhs[n+3][i][j][k];
      lhs[n+3][i][j1][k] = lhs[n+3][i][j1][k] -
	lhs[n+1][i][j1][k]*lhs[n+4][i][j][k];
      for (m = 0; m < 3; m++) {
	rhs[m][i][j1][k] = rhs[m][i][j1][k] -
	  lhs[n+1][i][j1][k]*rhs[m][i][j][k];
      }
      fac2               = 1./lhs[n+2][i][j1][k];
      for (m = 0; m < 3; m++) {
	rhs[m][i][j1][k] = fac2*rhs[m][i][j1][k];
      }
    }
  }
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5;
    for (j = 0; j <= grid_points[1]-3; j++) {
      j1 = j  + 1;
      j2 = j  + 2;
      for (i = 1; i <= grid_points[0]-2; i++) {
	for (k = 1; k <= grid_points[2]-2; k++) {
	  fac1               = 1./lhs[n+2][i][j][k];
	  lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k];
	  lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k];
	  rhs[m][i][j][k] = fac1*rhs[m][i][j][k];
	  lhs[n+2][i][j1][k] = lhs[n+2][i][j1][k] -
	    lhs[n+1][i][j1][k]*lhs[n+3][i][j][k];
	  lhs[n+3][i][j1][k] = lhs[n+3][i][j1][k] -
	    lhs[n+1][i][j1][k]*lhs[n+4][i][j][k];
	  rhs[m][i][j1][k] = rhs[m][i][j1][k] -
	    lhs[n+1][i][j1][k]*rhs[m][i][j][k];
	  lhs[n+1][i][j2][k] = lhs[n+1][i][j2][k] -
	    lhs[n+0][i][j2][k]*lhs[n+3][i][j][k];
	  lhs[n+2][i][j2][k] = lhs[n+2][i][j2][k] -
	    lhs[n+0][i][j2][k]*lhs[n+4][i][j][k];
	  rhs[m][i][j2][k] = rhs[m][i][j2][k] -
	    lhs[n+0][i][j2][k]*rhs[m][i][j][k];
	}
      }
    }
    j  = grid_points[1]-2;
    j1 = grid_points[1]-1;
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	fac1               = 1./lhs[n+2][i][j][k];
	lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k];
	lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k];
	rhs[m][i][j][k]     = fac1*rhs[m][i][j][k];
	lhs[n+2][i][j1][k] = lhs[n+2][i][j1][k] -
	  lhs[n+1][i][j1][k]*lhs[n+3][i][j][k];
	lhs[n+3][i][j1][k] = lhs[n+3][i][j1][k] -
	  lhs[n+1][i][j1][k]*lhs[n+4][i][j][k];
	rhs[m][i][j1][k]   = rhs[m][i][j1][k] -
	  lhs[n+1][i][j1][k]*rhs[m][i][j][k];
	fac2               = 1./lhs[n+2][i][j1][k];
	rhs[m][i][j1][k]   = fac2*rhs[m][i][j1][k];
      }
    }
  }
  j  = grid_points[1]-2;
  j1 = grid_points[1]-1;
  n = 0;
  for (m = 0; m < 3; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	rhs[m][i][j][k] = rhs[m][i][j][k] -
	  lhs[n+3][i][j][k]*rhs[m][i][j1][k];
      }
    }
  }
  for (m = 3; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	n = (m-3+1)*5;
	rhs[m][i][j][k] = rhs[m][i][j][k] -
	  lhs[n+3][i][j][k]*rhs[m][i][j1][k];
      }
    }
  }
  n = 0;
  for (m = 0; m < 3; m++) {
    for (j = grid_points[1]-3; j >= 0; j--) {
      j1 = j  + 1;
      j2 = j  + 2;
      for (i = 1; i <= grid_points[0]-2; i++) {
	for (k = 1; k <= grid_points[2]-2; k++) {
	  rhs[m][i][j][k] = rhs[m][i][j][k] - 
	    lhs[n+3][i][j][k]*rhs[m][i][j1][k] -
	    lhs[n+4][i][j][k]*rhs[m][i][j2][k];
	}
      }
    }
  }
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5;
    for (j = grid_points[1]-3; j >= 0; j--) {
      j1 = j  + 1;
      j2 = j1 + 1;
      for (i = 1; i <= grid_points[0]-2; i++) {
	for (k = 1; k <= grid_points[2]-2; k++) {
	  rhs[m][i][j][k] = rhs[m][i][j][k] - 
	    lhs[n+3][i][j][k]*rhs[m][i][j1][k] -
	    lhs[n+4][i][j][k]*rhs[m][i][j2][k];
	}
      }
    }
  }
}
  pinvr();
}