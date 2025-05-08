static void lhsz(void) {
  double ru1;
  int i, j, k;

  // The calculations for different (i, j) grid points are largely independent.
  // Grouping all k-dimension calculations for a single (i, j) pair together
  // makes the outer i and j loops ideal candidates for parallelization.
  // The temporary arrays 'cv' and 'rhos' are specific to each (i, j) pair.
  // Declaring them inside the j loop scope ensures they can be easily privatized
  // per thread when parallelizing the i or j loops.

  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {

      // Declare temporaries 'cv' and 'rhos' here. In an OpenMP parallel region
      // around this loop nest, these arrays would naturally be candidates
      // for privatization (e.g., using #pragma omp parallel for private(cv, rhos)).
      double cv[grid_points[2]];
      double rhos[grid_points[2]];

      // First k loop: Calculate cv and rhos for the current (i, j) slice.
      // This loop is independent across k iterations for fixed (i,j).
      for (k = 0; k <= grid_points[2]-1; k++) {
	ru1 = c3c4*rho_i[i][j][k];
	cv[k] = ws[i][j][k];
	rhos[k] = max(dz4 + con43 * ru1,
		      max(dz5 + c1c5 * ru1,
			  max(dzmax + ru1,
			      dz1)));
      }

      // Second k loop: Calculate initial lhs[0-4] for the current (i, j) slice.
      // This loop reads values from cv and rhos calculated in the previous loop
      // for the same (i, j). Writes are to lhs[band][i][j][k], which are
      // disjoint for different (i, j) pairs.
      for (k = 1; k <= grid_points[2]-2; k++) {
	lhs[0][i][j][k] =  0.0;
	lhs[1][i][j][k] = -dttz2 * cv[k-1] - dttz1 * rhos[k-1];
	lhs[2][i][j][k] =  1.0 + c2dttz1 * rhos[k];
	lhs[3][i][j][k] =  dttz2 * cv[k+1] - dttz1 * rhos[k+1];
	lhs[4][i][j][k] =  0.0;
      }

      // Boundary updates for k=1 and k=2 for the current (i, j) slice.
      // These modify lhs[0-4][i][j][k] at specific k values.
      // Writes are disjoint for different (i, j) pairs.
      k = 1;
      lhs[2][i][j][k] = lhs[2][i][j][k] + comz5;
      lhs[3][i][j][k] = lhs[3][i][j][k] - comz4;
      lhs[4][i][j][k] = lhs[4][i][j][k] + comz1;
      lhs[1][i][j][k+1] = lhs[1][i][j][k+1] - comz4;
      lhs[2][i][j][k+1] = lhs[2][i][j][k+1] + comz6;
      lhs[3][i][j][k+1] = lhs[3][i][j][k+1] - comz4;
      lhs[4][i][j][k+1] = lhs[4][i][j][k+1] + comz1;

      // Main k loop updates for k=3..grid_points[2]-4 for the current (i, j) slice.
      // This loop is independent across k iterations for fixed (i,j).
      // Writes are to lhs[band][i][j][k], disjoint for different (i,j).
      for (k = 3; k <= grid_points[2]-4; k++) {
	lhs[0][i][j][k] = lhs[0][i][j][k] + comz1;
	lhs[1][i][j][k] = lhs[1][i][j][k] - comz4;
	lhs[2][i][j][k] = lhs[2][i][j][k] + comz6;
	lhs[3][i][j][k] = lhs[3][i][j][k] - comz4;
	lhs[4][i][j][k] = lhs[4][i][j][k] + comz1;
      }

      // Boundary updates for k=grid_points[2]-3 and k=grid_points[2]-2
      // for the current (i, j) slice.
      // These modify lhs[0-4][i][j][k] at specific k values.
      // Writes are disjoint for different (i, j) pairs.
      k = grid_points[2]-3;
      lhs[0][i][j][k] = lhs[0][i][j][k] + comz1;
      lhs[1][i][j][k] = lhs[1][i][j][k] - comz4;
      lhs[2][i][j][k] = lhs[2][i][j][k] + comz6;
      lhs[3][i][j][k] = lhs[3][i][j][k] - comz4;
      lhs[0][i][j][k+1] = lhs[0][i][j][k+1] + comz1;
      lhs[1][i][j][k+1] = lhs[1][i][j][k+1] - comz4;
      lhs[2][i][j][k+1] = lhs[2][i][j][k+1] + comz5;

      // Final k loop: Calculate lhs[5-14] for the current (i, j) slice.
      // Reads the updated lhs[0-4][i][j][k] and speed[i][j][k].
      // Writes are to lhs[band][i][j][k], disjoint for different (i,j).
      for (k = 1; k <= grid_points[2]-2; k++) {
	lhs[0+5][i][j][k]  = lhs[0][i][j][k];
	lhs[1+5][i][j][k]  = lhs[1][i][j][k] -
	  dttz2 * speed[i][j][k-1];
	lhs[2+5][i][j][k]  = lhs[2][i][j][k];
	lhs[3+5][i][j][k]  = lhs[3][i][j][k] +
	  dttz2 * speed[i][j][k+1];
	lhs[4+5][i][j][k]  = lhs[4][i][j][k];
	lhs[0+10][i][j][k] = lhs[0][i][j][k];
	lhs[1+10][i][j][k] = lhs[1][i][j][k] +
	  dttz2 * speed[i][j][k-1];
	lhs[2+10][i][j][k] = lhs[2][i][j][k];
	lhs[3+10][i][j][k] = lhs[3][i][j][k] -
	  dttz2 * speed[i][j][k+1];
	lhs[4+10][i][j][k]  = lhs[4][i][j][k];
      }
    } // end j loop
  } // end i loop
}