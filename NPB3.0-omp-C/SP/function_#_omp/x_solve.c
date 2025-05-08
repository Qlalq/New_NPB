static void x_solve(void) {
{
  int i, j, k, n, i1, i2, m;
  double fac1, fac2;
  lhsx();
  n = 0;
  for (i = 0; i <= grid_points[0]-3; i++) {
    i1 = i  + 1;
    i2 = i  + 2;
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
  i  = grid_points[0]-2;
  i1 = grid_points[0]-1;
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
      fac2               = 1./lhs[n+2][i1][j][k];
      for (m = 0; m < 3; m++) {
	rhs[m][i1][j][k] = fac2*rhs[m][i1][j][k];
      }
    }
  }
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5;
    for (i = 0; i <= grid_points[0]-3; i++) {
      i1 = i  + 1;
      i2 = i  + 2;
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
    i  = grid_points[0]-2;
    i1 = grid_points[0]-1;
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
	fac2               = 1./lhs[n+2][i1][j][k];
	rhs[m][i1][j][k]   = fac2*rhs[m][i1][j][k];
      }
    }
  }
  i  = grid_points[0]-2;
  i1 = grid_points[0]-1;
  n = 0;
  for (m = 0; m < 3; m++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	rhs[m][i][j][k] = rhs[m][i][j][k] -
	  lhs[n+3][i][j][k]*rhs[m][i1][j][k];
      }
    }
  }
  for (m = 3; m < 5; m++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	n = (m-3+1)*5;
	rhs[m][i][j][k] = rhs[m][i][j][k] -
	  lhs[n+3][i][j][k]*rhs[m][i1][j][k];
      }
    }
  }
  n = 0;
  for (i = grid_points[0]-3; i >= 0; i--) {
    i1 = i  + 1;
    i2 = i  + 2;
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
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5;
    for (i = grid_points[0]-3; i >= 0; i--) {
      i1 = i  + 1;
      i2 = i  + 2;
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
  ninvr();
}