static void lhsinit(void) {

#pragma omp parallel
{
  int i, j, k, m, n;

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c     zero the whole left hand side for starters
c-------------------------------------------------------------------*/
#pragma omp for
  for (i = 0; i < grid_points[0]; i++) {
    for (j = 0; j < grid_points[1]; j++) {
      for (k = 0; k < grid_points[2]; k++) {
	for (m = 0; m < 5; m++) {
	  for (n = 0; n < 5; n++) {
	    lhs[i][j][k][0][m][n] = 0.0;
	    lhs[i][j][k][1][m][n] = 0.0;
	    lhs[i][j][k][2][m][n] = 0.0;
	  }
	}
      }
    }
  }

/*--------------------------------------------------------------------
c     next, set all diagonal values to 1. This is overkill, but convenient
c-------------------------------------------------------------------*/
#pragma omp for  
  for (i = 0; i < grid_points[0]; i++) {
    for (j = 0; j < grid_points[1]; j++) {
      for (k = 0; k < grid_points[2]; k++) {
	for (m = 0; m < 5; m++) {
	  lhs[i][j][k][1][m][m] = 1.0;
	}
      }
    }
  }
}

}