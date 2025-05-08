static void z_solve(void) {
  int i, j, k, n, k1, k2, m;
  double fac1, fac2;

  lhsz();

  // Block 1: Forward pass (m = 0 to 2), k = 0 to K-3
  // The loop over k has loop-carried dependencies. It must remain sequential.
  // The loops over i and j (and the inner loop over m) are independent for a fixed k.
  // These nested loops (i, j, m) can be parallelized.
  n = 0;
  for (k = 0; k <= grid_points[2]-3; k++) {
    k1 = k  + 1;
    k2 = k  + 2;
    // OpenMP Parallel region candidate: parallelize the i and j loops
    // #pragma omp parallel for collapse(2) private(i, j, fac1, m) if(...) // Add pragmas later
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
        // Calculations for (i, j, k) using values at k, k+1, k+2
        // Variables i, j, fac1, m are candidates for privatization
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

  // Block 2: Forward pass (m = 0 to 2), k = K-2, k1 = K-1 (Boundary)
  // These loops (i, j, m) are independent.
  k  = grid_points[2]-2;
  k1 = grid_points[2]-1;
  n = 0;
  // OpenMP Parallel region candidate: parallelize the i and j loops
  // #pragma omp parallel for collapse(2) private(i, j, fac1, fac2, m) if(...) // Add pragmas later
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      // Calculations for (i, j, k) using values at k, k1
      // Variables i, j, fac1, fac2, m are candidates for privatization
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

  // Block 3: Forward pass (m = 3 to 4), k = 0 to K-3
  // The loop over m is independent. It can be parallelized as the outermost loop.
  // Inside the m loop, the k loop is sequential due to dependency.
  // Inside the k loop, the i and j loops are independent and can be parallelized.
  // OpenMP Parallel region candidate: parallelize the m loop
  // #pragma omp parallel for private(m, n, k, k1, k2, fac1, i, j) if(...) // Add pragmas later
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5; // n depends on m, must be private or calculated inside m loop

    // Sequential loop over k due to data dependency
    for (k = 0; k <= grid_points[2]-3; k++) {
      k1 = k  + 1;
      k2 = k  + 2;
      // OpenMP Parallel region candidate: parallelize the i and j loops (within m and k)
      // #pragma omp parallel for collapse(2) private(i, j, fac1) if(...) // Add pragmas later
      for (i = 1; i <= grid_points[0]-2; i++) {
        for (j = 1; j <= grid_points[1]-2; j++) {
          // Calculations for (i, j, k) using values at k, k+1, k+2
          // Variables k, k1, k2, fac1, i, j are candidates for privatization
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
  }

  // Block 4: Forward pass (m = 3 to 4), k = K-2, k1 = K-1 (Boundary)
  // The loop over m is independent. It can be parallelized.
  // The loops over i and j are independent. They can be parallelized.
  k  = grid_points[2]-2;
  k1 = grid_points[2]-1;
  // OpenMP Parallel region candidate: parallelize the m loop, and/or the i/j loops
  // #pragma omp parallel for private(m, n, fac1, fac2, i, j) if(...) // Add pragmas later
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5; // n depends on m, must be private or calculated inside m loop

    // OpenMP Parallel region candidate: parallelize the i and j loops (within m)
    // #pragma omp parallel for collapse(2) private(i, j, fac1, fac2) if(...) // Add pragmas later
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
        // Calculations for (i, j, k) using values at k, k1
        // Variables fac1, fac2, i, j are candidates for privatization
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

  // Block 5: Backward pass (m = 0 to 2), k = K-2, k1 = K-1 (Boundary)
  // These loops (m, i, j) are independent. They can be parallelized.
  k  = grid_points[2]-2;
  k1 = grid_points[2]-1;
  n = 0;
  // OpenMP Parallel region candidate: parallelize the m, i, and j loops
  // #pragma omp parallel for collapse(3) private(m, i, j) if(...) // Add pragmas later
  for (m = 0; m < 3; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
        // Calculation for (i, j, k) using values at k and k1
        // Variables m, i, j are candidates for privatization
        rhs[m][i][j][k] = rhs[m][i][j][k] -
          lhs[n+3][i][j][k]*rhs[m][i][j][k1];
      }
    }
  }

  // Block 6: Backward pass (m = 3 to 4), k = K-2, k1 = K-1 (Boundary)
  // These loops (m, i, j) are independent. They can be parallelized.
  k  = grid_points[2]-2;
  k1 = grid_points[2]-1;
  // OpenMP Parallel region candidate: parallelize the m, i, and j loops
  // #pragma omp parallel for collapse(3) private(m, n, i, j) if(...) // Add pragmas later
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5; // n depends on m, must be private or calculated inside m loop
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
        // Calculation for (i, j, k) using values at k and k1
        // Variables m, n, i, j are candidates for privatization
        rhs[m][i][j][k] = rhs[m][i][j][k] -
          lhs[n+3][i][j][k]*rhs[m][i][j][k1];
      }
    }
  }

  // Block 7: Backward pass (m = 0 to 2), k = K-3 down to 0
  // The loop over k has loop-carried dependencies (counting down). It must remain sequential.
  // The loops over m, i, and j are independent for a fixed k.
  // These nested loops (m, i, j) can be parallelized.
  n = 0;
  for (k = grid_points[2]-3; k >= 0; k--) {
    k1 = k  + 1;
    k2 = k  + 2;
    // OpenMP Parallel region candidate: parallelize the m, i, and j loops
    // #pragma omp parallel for collapse(3) private(m, i, j, k1, k2) if(...) // Add pragmas later
    for (m = 0; m < 3; m++) {
      for (i = 1; i <= grid_points[0]-2; i++) {
        for (j = 1; j <= grid_points[1]-2; j++) {
          // Calculation for (i, j, k) using values at k, k1, k2
          // Variables m, i, j, k1, k2 are candidates for privatization
          rhs[m][i][j][k] = rhs[m][i][j][k] -
            lhs[n+3][i][j][k]*rhs[m][i][j][k1] -
            lhs[n+4][i][j][k]*rhs[m][i][j][k2];
        }
      }
    }
  }

  // Block 8: Backward pass (m = 3 to 4), k = K-3 down to 0
  // The loop over k has loop-carried dependencies (counting down). It must remain sequential.
  // The loops over m, i, and j are independent for a fixed k.
  // These nested loops (m, i, j) can be parallelized.
  for (k = grid_points[2]-3; k >= 0; k--) {
    k1 = k  + 1;
    k2 = k  + 2;
    // OpenMP Parallel region candidate: parallelize the m, i, and j loops
    // #pragma omp parallel for collapse(3) private(m, n, i, j, k1, k2) if(...) // Add pragmas later
    for (m = 3; m < 5; m++) {
      n = (m-3+1)*5; // n depends on m, must be private or calculated inside m loop
      for (i = 1; i <= grid_points[0]-2; i++) {
        for (j = 1; j <= grid_points[1]-2; j++) {
          // Calculation for (i, j, k) using values at k, k1, k2
          // Variables m, n, i, j, k1, k2 are candidates for privatization
          rhs[m][i][j][k] = rhs[m][i][j][k] -
            lhs[n+3][i][j][k]*rhs[m][i][j][k1] -
            lhs[n+4][i][j][k]*rhs[m][i][j][k2];
        }
      }
    }
  }

  tzetar();
}