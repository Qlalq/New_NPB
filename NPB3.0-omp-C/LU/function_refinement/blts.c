static void blts (int nx, int ny, int nz, int k,
		  double omega,
		  double v[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5],
		  double ldz[ISIZ1][ISIZ2][5][5],
		  double ldy[ISIZ1][ISIZ2][5][5],
		  double ldx[ISIZ1][ISIZ2][5][5],
		  double d[ISIZ1][ISIZ2][5][5],
		  int ist, int iend, int jst, int jend,
		  int nx0, int ny0 ) {
  int i, j, m;
  double tmp, tmp1;
  // tmat is local to the inner j loop, which is suitable for privatization in OpenMP
  double tmat[5][5]; 

  // First loop block: This loop calculates the contribution from the k-1 layer.
  // There are no loop-carried dependencies on i or j within this block.
  // This block is highly parallelizable across i and j.
  for (i = ist; i <= iend; i++) {
    for (j = jst; j <= jend; j++) {
      for (m = 0; m < 5; m++) {
	v[i][j][k][m] = v[i][j][k][m]
	  - omega * (  ldz[i][j][m][0] * v[i][j][k-1][0]
		       + ldz[i][j][m][1] * v[i][j][k-1][1]
		       + ldz[i][j][m][2] * v[i][j][k-1][2]
		       + ldz[i][j][m][3] * v[i][j][k-1][3]
		       + ldz[i][j][m][4] * v[i][j][k-1][4]  );
      }
    }
  }

  // Second loop block: This block solves the system involving dependencies
  // on i-1 and j-1 from the current k layer, and performs a local matrix solve.
  // The j loop has a loop-carried dependency (v[i][j] depends on v[i][j-1]).
  // The i loop has a dependency on i-1 (v[i][j] depends on v[i-1][j]).
  // This structure suggests a wavefront or pipeline parallelization pattern across 'i'.
  // The manual flag-based synchronization has been removed to prepare for
  // standard OpenMP dependency/synchronization primitives (e.g., using depend clauses).
  for (i = ist; i <= iend; i++) {
    // OpenMP synchronization/dependency directives would be placed here later
    // (e.g., #pragma omp task depend(in:v[i-1][jst:jend][k], ...) depend(out:v[i][jst:jend][k]))
    for (j = jst; j <= jend; j++) {
      // Apply dependencies from i-1 and j-1
      for (m = 0; m < 5; m++) {
	v[i][j][k][m] = v[i][j][k][m]
	  - omega * ( ldy[i][j][m][0] * v[i][j-1][k][0]
		      + ldx[i][j][m][0] * v[i-1][j][k][0]
		      + ldy[i][j][m][1] * v[i][j-1][k][1]
		      + ldx[i][j][m][1] * v[i-1][j][k][1]
		      + ldy[i][j][m][2] * v[i][j-1][k][2]
		      + ldx[i][j][m][2] * v[i-1][j][k][2]
		      + ldy[i][j][m][3] * v[i][j-1][k][3]
		      + ldx[i][j][m][3] * v[i-1][j][k][3]
		      + ldy[i][j][m][4] * v[i][j-1][k][4]
		      + ldx[i][j][m][4] * v[i-1][j][k][4] );
      }

      // Local matrix solve (Gaussian elimination and back substitution) for the current (i, j, k) cell.
      // This part is local and does not introduce new dependencies between (i,j) cells at the same k.
      for (m = 0; m < 5; m++) {
	tmat[m][0] = d[i][j][m][0];
	tmat[m][1] = d[i][j][m][1];
	tmat[m][2] = d[i][j][m][2];
	tmat[m][3] = d[i][j][m][3];
	tmat[m][4] = d[i][j][m][4];
      }
      tmp1 = 1.0 / tmat[0][0];
      tmp = tmp1 * tmat[1][0];
      tmat[1][1] =  tmat[1][1] - tmp * tmat[0][1];
      tmat[1][2] =  tmat[1][2] - tmp * tmat[0][2];
      tmat[1][3] =  tmat[1][3] - tmp * tmat[0][3];
      tmat[1][4] =  tmat[1][4] - tmp * tmat[0][4];
      v[i][j][k][1] = v[i][j][k][1] - v[i][j][k][0] * tmp;
      tmp = tmp1 * tmat[2][0];
      tmat[2][1] =  tmat[2][1] - tmp * tmat[0][1];
      tmat[2][2] =  tmat[2][2] - tmp * tmat[0][2];
      tmat[2][3] =  tmat[2][3] - tmp * tmat[0][3];
      tmat[2][4] =  tmat[2][4] - tmp * tmat[0][4];
      v[i][j][k][2] = v[i][j][k][2] - v[i][j][k][0] * tmp;
      tmp = tmp1 * tmat[3][0];
      tmat[3][1] =  tmat[3][1] - tmp * tmat[0][1];
      tmat[3][2] =  tmat[3][2] - tmp * tmat[0][2];
      tmat[3][3] =  tmat[3][3] - tmp * tmat[0][3];
      tmat[3][4] =  tmat[3][4] - tmp * tmat[0][4];
      v[i][j][k][3] = v[i][j][k][3] - v[i][j][k][0] * tmp;
      tmp = tmp1 * tmat[4][0];
      tmat[4][1] =  tmat[4][1] - tmp * tmat[0][1];
      tmat[4][2] =  tmat[4][2] - tmp * tmat[0][2];
      tmat[4][3] =  tmat[4][3] - tmp * tmat[0][3];
      tmat[4][4] =  tmat[4][4] - tmp * tmat[0][4];
      v[i][j][k][4] = v[i][j][k][4] - v[i][j][k][0] * tmp;
      tmp1 = 1.0 / tmat[ 1][1];
      tmp = tmp1 * tmat[ 2][1];
      tmat[2][2] =  tmat[2][2] - tmp * tmat[1][2];
      tmat[2][3] =  tmat[2][3] - tmp * tmat[1][3];
      tmat[2][4] =  tmat[2][4] - tmp * tmat[1][4];
      v[i][j][k][2] = v[i][j][k][2] - v[i][j][k][1] * tmp;
      tmp = tmp1 * tmat[3][1];
      tmat[3][2] =  tmat[3][2] - tmp * tmat[1][2];
      tmat[3][3] =  tmat[3][3] - tmp * tmat[1][3];
      tmat[3][4] =  tmat[3][4] - tmp * tmat[1][4];
      v[i][j][k][3] = v[i][j][k][3] - v[i][j][k][1] * tmp;
      tmp = tmp1 * tmat[4][1];
      tmat[4][2] =  tmat[4][2] - tmp * tmat[1][2];
      tmat[4][3] =  tmat[4][3] - tmp * tmat[1][3];
      tmat[4][4] =  tmat[4][4] - tmp * tmat[1][4];
      v[i][j][k][4] = v[i][j][k][4] - v[i][j][k][1] * tmp;
      tmp1 = 1.0 / tmat[2][2];
      tmp = tmp1 * tmat[3][2];
      tmat[3][3] =  tmat[3][3] - tmp * tmat[2][3];
      tmat[3][4] =  tmat[3][4] - tmp * tmat[2][4];
      v[i][j][k][3] = v[i][j][k][3]
        - v[i][j][k][2] * tmp;
      tmp = tmp1 * tmat[4][2];
      tmat[4][3] =  tmat[4][3] - tmp * tmat[2][3];
      tmat[4][4] =  tmat[4][4] - tmp * tmat[2][4];
      v[i][j][k][4] = v[i][j][k][4]
	- v[i][j][k][2] * tmp;
      tmp1 = 1.0 / tmat[3][3];
      tmp = tmp1 * tmat[4][3];
      tmat[4][4] =  tmat[4][4] - tmp * tmat[3][4];
      v[i][j][k][4] = v[i][j][k][4]
	- v[i][j][k][3] * tmp;
      v[i][j][k][4] = v[i][j][k][4]
	/ tmat[4][4];
      v[i][j][k][3] = v[i][j][k][3]
	- tmat[3][4] * v[i][j][k][4];
      v[i][j][k][3] = v[i][j][k][3]
	/ tmat[3][3];
      v[i][j][k][2] = v[i][j][k][2]
	- tmat[2][3] * v[i][j][k][3]
	- tmat[2][4] * v[i][j][k][4];
      v[i][j][k][2] = v[i][j][k][2]
	/ tmat[2][2];
      v[i][j][k][1] = v[i][j][k][1]
	- tmat[1][2] * v[i][j][k][2]
	- tmat[1][3] * v[i][j][k][3]
	- tmat[1][4] * v[i][j][k][4];
      v[i][j][k][1] = v[i][j][k][1]
	/ tmat[1][1];
      v[i][j][k][0] = v[i][j][k][0]
	- tmat[0][1] * v[i][j][k][1]
	- tmat[0][2] * v[i][j][k][2]
	- tmat[0][3] * v[i][j][k][3]
	- tmat[0][4] * v[i][j][k][4];
      v[i][j][k][0] = v[i][j][k][0]
	/ tmat[0][0];
    }
    // OpenMP signaling logic would be placed here later
    // (e.g., #pragma omp task depend(out:v[i][jst:jend][k]))
  }
}