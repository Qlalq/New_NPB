static void buts(int nx, int ny, int nz, int k,
		 double omega,
/*--------------------------------------------------------------------
c   To improve cache performance, second two dimensions padded by 1 
c   for even number sizes only.  Only needed in v.
--------------------------------------------------------------------*/
		 double v[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5],
		 double tv[ISIZ1][ISIZ2][5],
		 double d[ISIZ1][ISIZ2][5][5],
		 double udx[ISIZ1][ISIZ2][5][5],
		 double udy[ISIZ1][ISIZ2][5][5],
		 double udz[ISIZ1][ISIZ2][5][5],
		 int ist, int iend, int jst, int jend,
		 int nx0, int ny0 ) {
/*--------------------------------------------------------------------
c
c   compute the regular-sparse, block upper triangular solution:
c
c                     v <-- ( U-inv ) * v
c
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j, m;
  double tmp, tmp1;
  double tmat[5][5];

#pragma omp for nowait schedule(static)
  for (i = iend; i >= ist; i--) {
    for (j = jend; j >= jst; j--) {
      for (m = 0; m < 5; m++) {
	tv[i][j][m] = 
	  omega * (  udz[i][j][m][0] * v[i][j][k+1][0]
		     + udz[i][j][m][1] * v[i][j][k+1][1]
		     + udz[i][j][m][2] * v[i][j][k+1][2]
		     + udz[i][j][m][3] * v[i][j][k+1][3]
		     + udz[i][j][m][4] * v[i][j][k+1][4] );
      }
    }
  }

#pragma omp for nowait schedule(static)
  for (i = iend; i >= ist; i--) {
#if defined(_OPENMP)      
    if (i != iend) {
      while (flag[i+1] == 0) {
#pragma omp flush(flag)
	;
      }
    }
    if (i != ist) {
      while (flag[i] == 1) {
#pragma omp flush(flag)
	;
      }
    }
#endif /* _OPENMP */
    
    for (j = jend; j >= jst; j--) {
      for (m = 0; m < 5; m++) {
	tv[i][j][m] = tv[i][j][m]
	  + omega * ( udy[i][j][m][0] * v[i][j+1][k][0]
		      + udx[i][j][m][0] * v[i+1][j][k][0]
		      + udy[i][j][m][1] * v[i][j+1][k][1]
		      + udx[i][j][m][1] * v[i+1][j][k][1]
		      + udy[i][j][m][2] * v[i][j+1][k][2]
		      + udx[i][j][m][2] * v[i+1][j][k][2]
		      + udy[i][j][m][3] * v[i][j+1][k][3]
		      + udx[i][j][m][3] * v[i+1][j][k][3]
		      + udy[i][j][m][4] * v[i][j+1][k][4]
		      + udx[i][j][m][4] * v[i+1][j][k][4] );
      }

/*--------------------------------------------------------------------
c   diagonal block inversion
--------------------------------------------------------------------*/
      for (m = 0; m < 5; m++) {
	tmat[m][0] = d[i][j][m][0];
	tmat[m][1] = d[i][j][m][1];
	tmat[m][2] = d[i][j][m][2];
	tmat[m][3] = d[i][j][m][3];
	tmat[m][4] = d[i][j][m][4];
      }

      tmp1 = 1.0 / tmat[0][0];
      tmp = tmp1 * tmat[1][0];
      tmat[1][1] =  tmat[1][1]
	- tmp * tmat[0][1];
      tmat[1][2] =  tmat[1][2]
	- tmp * tmat[0][2];
      tmat[1][3] =  tmat[1][3]
	- tmp * tmat[0][3];
      tmat[1][4] =  tmat[1][4]
	- tmp * tmat[0][4];
      tv[i][j][1] = tv[i][j][1]
	- tv[i][j][0] * tmp;

      tmp = tmp1 * tmat[2][0];
      tmat[2][1] =  tmat[2][1]
	- tmp * tmat[0][1];
      tmat[2][2] =  tmat[2][2]
	- tmp * tmat[0][2];
      tmat[2][3] =  tmat[2][3]
	- tmp * tmat[0][3];
      tmat[2][4] =  tmat[2][4]
	- tmp * tmat[0][4];
      tv[i][j][2] = tv[i][j][2]
	- tv[i][j][0] * tmp;

      tmp = tmp1 * tmat[3][0];
      tmat[3][1] =  tmat[3][1]
	- tmp * tmat[0][1];
      tmat[3][2] =  tmat[3][2]
	- tmp * tmat[0][2];
      tmat[3][3] =  tmat[3][3]
	- tmp * tmat[0][3];
      tmat[3][4] =  tmat[3][4]
	- tmp * tmat[0][4];
      tv[i][j][3] = tv[i][j][3]
	- tv[i][j][0] * tmp;

      tmp = tmp1 * tmat[4][0];
      tmat[4][1] =  tmat[4][1]
	- tmp * tmat[0][1];
      tmat[4][2] =  tmat[4][2]
	- tmp * tmat[0][2];
      tmat[4][3] =  tmat[4][3]
	- tmp * tmat[0][3];
      tmat[4][4] =  tmat[4][4]
	- tmp * tmat[0][4];
      tv[i][j][4] = tv[i][j][4]
	- tv[i][j][0] * tmp;

      tmp1 = 1.0 / tmat[1][1];
      tmp = tmp1 * tmat[2][1];
      tmat[2][2] =  tmat[2][2]
	- tmp * tmat[1][2];
      tmat[2][3] =  tmat[2][3]
	- tmp * tmat[1][3];
      tmat[2][4] =  tmat[2][4]
	- tmp * tmat[1][4];
      tv[i][j][2] = tv[i][j][2]
	- tv[i][j][1] * tmp;

      tmp = tmp1 * tmat[3][1];
      tmat[3][2] =  tmat[3][2]
	- tmp * tmat[1][2];
      tmat[3][3] =  tmat[3][3]
	- tmp * tmat[1][3];
      tmat[3][4] =  tmat[3][4]
	- tmp * tmat[1][4];
      tv[i][j][3] = tv[i][j][3]
	- tv[i][j][1] * tmp;

      tmp = tmp1 * tmat[4][1];
      tmat[4][2] =  tmat[4][2]
	- tmp * tmat[1][2];
      tmat[4][3] =  tmat[4][3]
	- tmp * tmat[1][3];
      tmat[4][4] =  tmat[4][4]
	- tmp * tmat[1][4];
      tv[i][j][4] = tv[i][j][4]
	- tv[i][j][1] * tmp;

      tmp1 = 1.0 / tmat[2][2];
      tmp = tmp1 * tmat[3][2];
      tmat[3][3] =  tmat[3][3]
	- tmp * tmat[2][3];
      tmat[3][4] =  tmat[3][4]
	- tmp * tmat[2][4];
      tv[i][j][3] = tv[i][j][3]
	- tv[i][j][2] * tmp;

      tmp = tmp1 * tmat[4][2];
      tmat[4][3] =  tmat[4][3]
	- tmp * tmat[2][3];
      tmat[4][4] =  tmat[4][4]
	- tmp * tmat[2][4];
      tv[i][j][4] = tv[i][j][4]
	- tv[i][j][2] * tmp;

      tmp1 = 1.0 / tmat[3][3];
      tmp = tmp1 * tmat[4][3];
      tmat[4][4] =  tmat[4][4]
	- tmp * tmat[3][4];
      tv[i][j][4] = tv[i][j][4]
	- tv[i][j][3] * tmp;

/*--------------------------------------------------------------------
c   back substitution
--------------------------------------------------------------------*/
      tv[i][j][4] = tv[i][j][4]
	/ tmat[4][4];

      tv[i][j][3] = tv[i][j][3]
	- tmat[3][4] * tv[i][j][4];
      tv[i][j][3] = tv[i][j][3]
	/ tmat[3][3];

      tv[i][j][2] = tv[i][j][2]
	- tmat[2][3] * tv[i][j][3]
	- tmat[2][4] * tv[i][j][4];
      tv[i][j][2] = tv[i][j][2]
	/ tmat[2][2];

      tv[i][j][1] = tv[i][j][1]
	- tmat[1][2] * tv[i][j][2]
	- tmat[1][3] * tv[i][j][3]
	- tmat[1][4] * tv[i][j][4];
      tv[i][j][1] = tv[i][j][1]
	/ tmat[1][1];

      tv[i][j][0] = tv[i][j][0]
	- tmat[0][1] * tv[i][j][1]
	- tmat[0][2] * tv[i][j][2]
	- tmat[0][3] * tv[i][j][3]
	- tmat[0][4] * tv[i][j][4];
      tv[i][j][0] = tv[i][j][0]
	/ tmat[0][0];

      v[i][j][k][0] = v[i][j][k][0] - tv[i][j][0];
      v[i][j][k][1] = v[i][j][k][1] - tv[i][j][1];
      v[i][j][k][2] = v[i][j][k][2] - tv[i][j][2];
      v[i][j][k][3] = v[i][j][k][3] - tv[i][j][3];
      v[i][j][k][4] = v[i][j][k][4] - tv[i][j][4];
    }
    
#if defined(_OPENMP)    
    if (i != iend) flag[i+1] = 0;
    if (i != ist) flag[i] = 1;
#pragma omp flush(flag)
#endif /* _OPENMP */    
  }
}