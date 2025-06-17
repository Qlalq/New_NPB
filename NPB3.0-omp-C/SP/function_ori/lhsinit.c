static void lhsinit(void) {

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

  int i, j, k, n;

/*--------------------------------------------------------------------
c     zap the whole left hand side for starters
c-------------------------------------------------------------------*/
  for (n = 0; n < 15; n++) {
#pragma omp for nowait
    for (i = 0; i < grid_points[0]; i++) {
      for (j = 0; j < grid_points[1]; j++) {
	for (k = 0; k < grid_points[2]; k++) {
	  lhs[n][i][j][k] = 0.0;
	}
      }
    }
  }
#pragma omp barrier  

/*--------------------------------------------------------------------
c      next, set all diagonal values to 1. This is overkill, but 
c      convenient
c-------------------------------------------------------------------*/
  for (n = 0; n < 3; n++) {
#pragma omp for    
    for (i = 0; i < grid_points[0]; i++) {
      for (j = 0; j < grid_points[1]; j++) {
	for (k = 0; k < grid_points[2]; k++) {
	  lhs[5*n+2][i][j][k] = 1.0;
	}
      }
    }
  }
}