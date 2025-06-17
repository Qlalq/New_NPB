static void ninvr(void) {

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c   block-diagonal matrix-vector multiplication              
c-------------------------------------------------------------------*/

  int i, j, k;
  double r1, r2, r3, r4, r5, t1, t2;
#pragma omp parallel for default(shared) private(i,j,k,r1,r2,r3,r4,r5,t1,t2)
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {

	r1 = rhs[0][i][j][k];
	r2 = rhs[1][i][j][k];
	r3 = rhs[2][i][j][k];
	r4 = rhs[3][i][j][k];
	r5 = rhs[4][i][j][k];
               
	t1 = bt * r3;
	t2 = 0.5 * ( r4 + r5 );

	rhs[0][i][j][k] = -r2;
	rhs[1][i][j][k] =  r1;
	rhs[2][i][j][k] = bt * ( r4 - r5 );
	rhs[3][i][j][k] = -t1 + t2;
	rhs[4][i][j][k] =  t1 + t2;
      }
    }
  }
}