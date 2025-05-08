static void z_solve(void) {
{
  int i, j, k, n, k1, k2, m;
  double fac1, fac2;
  lhsz();
  n = 0;
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
  k  = grid_points[2]-2;
  k1 = grid_points[2]-1;
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
      fac2               = 1./lhs[n+2][i][j][k1];
      for (m = 0; m < 3; m++) {
	rhs[m][i][j][k1] = fac2*rhs[m][i][j][k1];
      }
    }
  }
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5;
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
    k  = grid_points[2]-2;
    k1 = grid_points[2]-1;
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
	fac2               = 1./lhs[n+2][i][j][k1];
	rhs[m][i][j][k1]   = fac2*rhs[m][i][j][k1];
      }
    }
  }
  k  = grid_points[2]-2;
  k1 = grid_points[2]-1;
  n = 0;
  for (m = 0; m < 3; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	rhs[m][i][j][k] = rhs[m][i][j][k] -
	  lhs[n+3][i][j][k]*rhs[m][i][j][k1];
      }
    }
  }
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5;
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	rhs[m][i][j][k] = rhs[m][i][j][k] -
	  lhs[n+3][i][j][k]*rhs[m][i][j][k1];
      }
    }
  }
  n = 0;
  for (m = 0; m < 3; m++) {
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
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5;
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