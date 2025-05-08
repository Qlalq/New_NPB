static void buts_refined(int nx, int ny, int nz, int k,
		 double omega,
		 double v[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5],
		 double tv[ISIZ1][ISIZ2][5],
		 double d[ISIZ1][ISIZ2][5][5],
		 double udx[ISIZ1][ISIZ2][5][5],
		 double udy[ISIZ1][ISIZ2][5][5],
		 double udz[ISIZ1][ISIZ2][5][5],
		 int ist, int iend, int jst, int jend,
		 int nx0, int ny0 ) {
  int i, j, m;

  // Part 1: Compute initial tv component (independent loop)
  // This loop iterates backward over i and j, calculating a part of tv
  // based on v at plane k+1. It does not write to v and its calculations
  // for each (i,j) pair are independent of other (i',j') pairs within this loop.
  // This loop is highly parallelizable over i and/or j.
  for (i = iend; i >= ist; i--) { // Independent i loop
    for (j = jend; j >= jst; j--) { // Independent j loop
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

  // Part 2: Dependent sweep (Backward i, Backward j)
  // This loop updates tv based on neighbors v[i+1][j][k] and v[i][j+1][k],
  // then performs a local 5x5 linear system solve (Gaussian elimination + back substitution)
  // using d[i][j] and the updated tv[i][j], and finally updates v[i][j][k].
  // The i loop has a loop-carried dependency because the calculation for (i,j)
  // depends on the value of v updated by iteration (i+1, j).
  // The j loop also has a loop-carried dependency because the calculation for (i,j)
  // depends on the value of v updated by iteration (i, j+1).
  // Direct parallelization of the i or j loop requires careful synchronization
  // (e.g., wavefront parallelization) or using OpenMP task dependencies.
  // The original non-standard OpenMP flag synchronization is removed here.
  for (i = iend; i >= ist; i--) { // Dependent i loop (on i+1)
    for (j = jend; j >= jst; j--) { // Dependent j loop (on j+1)

      // Update tv[i][j] using neighbors v[i+1][j][k] and v[i][j+1][k]
      // This calculation reads from v at i+1 and j+1 for the current k plane.
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

      // Perform 5x5 Gaussian elimination and back substitution for (i,j)
      // This part solves a local linear system for delta_v[i][j] using the matrix d[i][j]
      // and the current tv[i][j] as the right-hand side.
      // This computation block for a given (i,j) pair is independent of other (i',j')
      // pairs once the input tv[i][j] is ready.
      // The variables tmat_local, tv_local, tmp, and tmp1 are local to this scope
      // and can be privatized in parallel execution.
      double tmat_local[5][5]; // Local copy for the matrix
      double tv_local[5];      // Local copy for the right-hand side vector
      double tmp, tmp1;        // Local temporary variables

      // Copy data to local variables for the solve
      for (m = 0; m < 5; m++) {
        tmat_local[m][0] = d[i][j][m][0];
        tmat_local[m][1] = d[i][j][m][1];
        tmat_local[m][2] = d[i][j][m][2];
        tmat_local[m][3] = d[i][j][m][3];
        tmat_local[m][4] = d[i][j][m][4];
        tv_local[m] = tv[i][j][m]; // Use the tv updated with neighbor terms
      }

      // Gaussian elimination (Forward elimination) on tmat_local and tv_local
      tmp1 = 1.0 / tmat_local[0][0];
      tmp = tmp1 * tmat_local[1][0];
      tmat_local[1][1] =  tmat_local[1][1] - tmp * tmat_local[0][1];
      tmat_local[1][2] =  tmat_local[1][2] - tmp * tmat_local[0][2];
      tmat_local[1][3] =  tmat_local[1][3] - tmp * tmat_local[0][3];
      tmat_local[1][4] =  tmat_local[1][4] - tmp * tmat_local[0][4];
      tv_local[1] = tv_local[1] - tv_local[0] * tmp;
      tmp = tmp1 * tmat_local[2][0];
      tmat_local[2][1] =  tmat_local[2][1] - tmp * tmat_local[0][1];
      tmat_local[2][2] =  tmat_local[2][2] - tmp * tmat_local[0][2];
      tmat_local[2][3] =  tmat_local[2][3] - tmp * tmat_local[0][3];
      tmat_local[2][4] =  tmat_local[2][4] - tmp * tmat_local[0][4];
      tv_local[2] = tv_local[2] - tv_local[0] * tmp;
      tmp = tmp1 * tmat_local[3][0];
      tmat_local[3][1] =  tmat_local[3][1] - tmp * tmat_local[0][1];
      tmat_local[3][2] =  tmat_local[3][2] - tmp * tmat_local[0][2];
      tmat_local[3][3] =  tmat_local[3][3] - tmp * tmat_local[0][3];
      tmat_local[3][4] =  tmat_local[3][4] - tmp * tmat_local[0][4];
      tv_local[3] = tv_local[3] - tv_local[0] * tmp;
      tmp = tmp1 * tmat_local[4][0];
      tmat_local[4][1] =  tmat_local[4][1] - tmp * tmat_local[0][1];
      tmat_local[4][2] =  tmat_local[4][2] - tmp * tmat_local[0][2];
      tmat_local[4][3] =  tmat_local[4][3] - tmp * tmat_local[0][3];
      tmat_local[4][4] =  tmat_local[4][4] - tmp * tmat_local[0][4];
      tv_local[4] = tv_local[4] - tv_local[0] * tmp;
      tmp1 = 1.0 / tmat_local[1][1];
      tmp = tmp1 * tmat_local[2][1];
      tmat_local[2][2] =  tmat_local[2][2] - tmp * tmat_local[1][2];
      tmat_local[2][3] =  tmat_local[2][3] - tmp * tmat_local[1][3];
      tmat_local[2][4] =  tmat_local[2][4] - tmp * tmat_local[1][4];
      tv_local[2] = tv_local[2] - tv_local[1] * tmp;
      tmp = tmp1 * tmat_local[3][1];
      tmat_local[3][2] =  tmat_local[3][2] - tmp * tmat_local[1][2];
      tmat_local[3][3] =  tmat_local[3][3] - tmp * tmat_local[1][3];
      tmat_local[3][4] =  tmat_local[3][4] - tmp * tmat_local[1][4];
      tv_local[3] = tv_local[3] - tv_local[1] * tmp;
      tmp = tmp1 * tmat_local[4][1];
      tmat_local[4][2] =  tmat_local[4][2] - tmp * tmat_local[1][2];
      tmat_local[4][3] =  tmat_local[4][3] - tmp * tmat_local[1][3];
      tmat_local[4][4] =  tmat_local[4][4] - tmp * tmat_local[1][4];
      tv_local[4] = tv_local[4] - tv_local[1] * tmp;
      tmp1 = 1.0 / tmat_local[2][2];
      tmp = tmp1 * tmat_local[3][2];
      tmat_local[3][3] =  tmat_local[3][3] - tmp * tmat_local[2][3];
      tmat_local[3][4] =  tmat_local[3][4] - tmp * tmat_local[2][4];
      tv_local[3] = tv_local[3] - tv_local[2] * tmp;
      tmp = tmp1 * tmat_local[4][2];
      tmat_local[4][3] =  tmat_local[4][3] - tmp * tmat_local[2][3];
      tmat_local[4][4] =  tmat_local[4][4] - tmp * tmat_local[2][4];
      tv_local[4] = tv_local[4] - tv_local[2] * tmp;
      tmp1 = 1.0 / tmat_local[3][3];
      tmp = tmp1 * tmat_local[4][3];
      tmat_local[4][4] =  tmat_local[4][4] - tmp * tmat_local[3][4];
      tv_local[4] = tv_local[4] - tv_local[3] * tmp;

      // Back substitution
      tv_local[4] = tv_local[4] / tmat_local[4][4];
      tv_local[3] = tv_local[3]
	- tmat_local[3][4] * tv_local[4];
      tv_local[3] = tv_local[3] / tmat_local[3][3];
      tv_local[2] = tv_local[2]
	- tmat_local[2][3] * tv_local[3]
	- tmat_local[2][4] * tv_local[4];
      tv_local[2] = tv_local[2] / tmat_local[2][2];
      tv_local[1] = tv_local[1]
	- tmat_local[1][2] * tv_local[2]
	- tmat_local[1][3] * tv_local[3]
	- tmat_local[1][4] * tv_local[4];
      tv_local[1] = tv_local[1] / tmat_local[1][1];
      tv_local[0] = tv_local[0]
	- tmat_local[0][1] * tv_local[1]
	- tmat_local[0][2] * tv_local[2]
	- tmat_local[0][3] * tv_local[3]
	- tmat_local[0][4] * tv_local[4];
      tv_local[0] = tv_local[0] / tmat_local[0][0];


      // Update v[i][j][k] using the solved tv_local (which represents delta_v)
      // This writes to v[i][j][k]. The value written here is read by iterations
      // (i-1, j) and (i, j-1) in the next loop iterations, creating the dependency.
      v[i][j][k][0] = v[i][j][k][0] - tv_local[0];
      v[i][j][k][1] = v[i][j][k][1] - tv_local[1];
      v[i][j][k][2] = v[i][j][k][2] - tv_local[2];
      v[i][j][k][3] = v[i][j][k][3] - tv_local[3];
      v[i][j][k][4] = v[i][j][k][4] - tv_local[4];

    } // end j loop
  } // end i loop
}