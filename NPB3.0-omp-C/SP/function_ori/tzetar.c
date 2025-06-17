static void tzetar(void) {

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c   block-diagonal matrix-vector multiplication                       
c-------------------------------------------------------------------*/

  int i, j, k;
  double t1, t2, t3, ac, xvel, yvel, zvel, r1, r2, r3, 
    r4, r5, btuz, acinv, ac2u, uzik1;
  
#pragma omp for private(i,j,k,t1,t2,t3,ac,xvel,yvel,zvel,r1,r2,r3,r4,r5,btuz,ac2u,uzik1)
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {

	xvel = us[i][j][k];
	yvel = vs[i][j][k];
	zvel = ws[i][j][k];
	ac   = speed[i][j][k];
	acinv = ainv[i][j][k];

	ac2u = ac*ac;

	r1 = rhs[0][i][j][k];
	r2 = rhs[1][i][j][k];
	r3 = rhs[2][i][j][k];
	r4 = rhs[3][i][j][k];
	r5 = rhs[4][i][j][k];

	uzik1 = u[0][i][j][k];
	btuz  = bt * uzik1;

	t1 = btuz*acinv * (r4 + r5);
	t2 = r3 + t1;
	t3 = btuz * (r4 - r5);

	rhs[0][i][j][k] = t2;
	rhs[1][i][j][k] = -uzik1*r2 + xvel*t2;
	rhs[2][i][j][k] =  uzik1*r1 + yvel*t2;
	rhs[3][i][j][k] =  zvel*t2  + t3;
	rhs[4][i][j][k] =  uzik1*(-xvel*r2 + yvel*r1) + 
	  qs[i][j][k]*t2 + c2iv*ac2u*t1 + zvel*t3;
      }
    }
  }
}