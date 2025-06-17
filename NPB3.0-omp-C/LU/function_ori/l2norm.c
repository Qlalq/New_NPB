static void l2norm (int nx0, int ny0, int nz0,
		    int ist, int iend,
		    int jst, int jend,
/*--------------------------------------------------------------------
c   To improve cache performance, second two dimensions padded by 1 
c   for even number sizes only.  Only needed in v.
--------------------------------------------------------------------*/
		    double v[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5],
		    double sum[5]) {

#pragma omp parallel 
{

/*--------------------------------------------------------------------
c   to compute the l2-norm of vector v.
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j, k, m;
  double sum0=0.0, sum1=0.0, sum2=0.0, sum3=0.0, sum4=0.0;

#pragma omp single  
  for (m = 0; m < 5; m++) {
    sum[m] = 0.0;
  }

#pragma omp for nowait
  for (i = ist; i <= iend; i++) {
    for (j = jst; j <= jend; j++) {
      for (k = 1; k <= nz0-2; k++) {
	  sum0 = sum0 + v[i][j][k][0] * v[i][j][k][0];
	  sum1 = sum1 + v[i][j][k][1] * v[i][j][k][1];
	  sum2 = sum2 + v[i][j][k][2] * v[i][j][k][2];
	  sum3 = sum3 + v[i][j][k][3] * v[i][j][k][3];
	  sum4 = sum4 + v[i][j][k][4] * v[i][j][k][4];
      }
    }
  }

#pragma omp critical
  {
      sum[0] += sum0;
      sum[1] += sum1;
      sum[2] += sum2;
      sum[3] += sum3;
      sum[4] += sum4;
  }
#pragma omp barrier  
  
#pragma omp single  
  for (m = 0;  m < 5; m++) {
    sum[m] = sqrt ( sum[m] / ( (nx0-2)*(ny0-2)*(nz0-2) ) );
  }
}
}