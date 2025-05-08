static void y_solve(void) {
  int i, j, k, n, j1, j2, m;
  double fac1, fac2;

  //-------------------------------------------------------------------------
  // Phase 1: Initialize LHS (implicitly done by lhsy) and
  //          perform forward elimination step along the y-direction (j).
  //          This process updates lhs and rhs arrays.
  //          The outer j-loop is sequential due to data dependencies.
  //          The inner i and k loops operate on independent grid points
  //          within a given set of j-planes and are parallelizable.
  //-------------------------------------------------------------------------
  lhsy();

  // Forward elimination for m = 0, 1, 2
  n = 0; // n is fixed at 0 for m = 0, 1, 2
  // Sequential loop over j due to dependencies (updates at j+1, j+2 depend on j)
  for (j = 0; j <= grid_points[1]-3; j++) {
    j1 = j  + 1;
    j2 = j  + 2;
    // These nested loops iterate over the i-k plane.
    // The calculations for different (i, k) pairs are independent.
    // This is a prime target for OpenMP parallelization (e.g., #pragma omp for collapse(2))
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
        // Calculations for m = 0, 1, 2 are inside the i-k loops
        fac1               = 1./lhs[n+2][i][j][k]; // uses j
        lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k]; // updates j
        lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k]; // updates j
        for (m = 0; m < 3; m++) {
          rhs[m][i][j][k] = fac1*rhs[m][i][j][k]; // updates j
        }
        // Updates at j+1 (j1) use previously updated values at j
        lhs[n+2][i][j1][k] = lhs[n+2][i][j1][k] -
          lhs[n+1][i][j1][k]*lhs[n+3][i][j][k];
        lhs[n+3][i][j1][k] = lhs[n+3][i][j1][k] -
          lhs[n+1][i][j1][k]*lhs[n+4][i][j][k];
        for (m = 0; m < 3; m++) {
          rhs[m][i][j1][k] = rhs[m][i][j1][k] -
            lhs[n+1][i][j1][k]*rhs[m][i][j][k];
        }
        // Updates at j+2 (j2) use previously updated values at j
        lhs[n+1][i][j2][k] = lhs[n+1][i][j2][k] -
          lhs[n+0][i][j2][k]*lhs[n+3][i][j][k];
        lhs[n+2][i][j2][k] = lhs[n+2][i][j2][k] -
          lhs[n+0][i][j2][k]*lhs[n+4][i][j][k];
        for (m = 0; m < 3; m++) {
          rhs[m][i][j2][k] = rhs[m][i][j2][k] -
            lhs[n+0][i][j2][k]*rhs[m][i][j][k];
        }
      }
    }
  }

  // Forward elimination boundary case for j = grid_points[1]-2 and j1 = grid_points[1]-1 (m = 0, 1, 2)
  j  = grid_points[1]-2;
  j1 = grid_points[1]-1;
  // These nested loops iterate over the i-k plane. Calculations for different (i, k) are independent.
  // This is also a target for OpenMP parallelization.
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (k = 1; k <= grid_points[2]-2; k++) {
      // Calculations for m = 0, 1, 2 are inside the i-k loops
      fac1               = 1./lhs[n+2][i][j][k]; // uses j
      lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k]; // updates j
      lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k]; // updates j
      for (m = 0; m < 3; m++) {
        rhs[m][i][j][k] = fac1*rhs[m][i][j][k]; // updates j
      }
      // Updates at j+1 (j1) use previously updated values at j
      lhs[n+2][i][j1][k] = lhs[n+2][i][j1][k] -
        lhs[n+1][i][j1][k]*lhs[n+3][i][j][k];
      lhs[n+3][i][j1][k] = lhs[n+3][i][j1][k] -
        lhs[n+1][i][j1][k]*lhs[n+4][i][j][k];
      for (m = 0; m < 3; m++) {
        rhs[m][i][j1][k] = rhs[m][i][j1][k] -
          lhs[n+1][i][j1][k]*rhs[m][i][j][k];
      }
      // Special handling for the last plane (j1) after update
      fac2               = 1./lhs[n+2][i][j1][k]; // uses j1
      for (m = 0; m < 3; m++) {
        rhs[m][i][j1][k] = fac2*rhs[m][i][j1][k]; // updates j1
      }
    }
  }

  // Forward elimination for m = 3, 4
  // Sequential loop over m
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5; // n depends on m

    // Sequential loop over j due to dependencies (updates at j+1, j+2 depend on j)
    for (j = 0; j <= grid_points[1]-3; j++) {
      j1 = j  + 1;
      j2 = j  + 2;
      // These nested loops iterate over the i-k plane.
      // Calculations for different (i, k) pairs are independent for a fixed m and j.
      // This is a prime target for OpenMP parallelization (e.g., #pragma omp for collapse(2))
      for (i = 1; i <= grid_points[0]-2; i++) {
        for (k = 1; k <= grid_points[2]-2; k++) {
          // Calculations for m = 3, 4 are inside the i-k loops
          fac1               = 1./lhs[n+2][i][j][k]; // uses j, n (depends on m)
          lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k]; // updates j, n
          lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k]; // updates j, n
          rhs[m][i][j][k] = fac1*rhs[m][i][j][k]; // updates j, m
          // Updates at j+1 (j1) use previously updated values at j
          lhs[n+2][i][j1][k] = lhs[n+2][i][j1][k] -
            lhs[n+1][i][j1][k]*lhs[n+3][i][j][k]; // uses j1, n, updated j, n
          lhs[n+3][i][j1][k] = lhs[n+3][i][j1][k] -
            lhs[n+1][i][j1][k]*lhs[n+4][i][j][k]; // uses j1, n, updated j, n
          rhs[m][i][j1][k] = rhs[m][i][j1][k] -
            lhs[n+1][i][j1][k]*rhs[m][i][j][k]; // uses j1, n, updated j, m
          // Updates at j+2 (j2) use previously updated values at j
          lhs[n+1][i][j2][k] = lhs[n+1][i][j2][k] -
            lhs[n+0][i][j2][k]*lhs[n+3][i][j][k]; // uses j2, n, updated j, n
          lhs[n+2][i][j2][k] = lhs[n+2][i][j2][k] -
            lhs[n+0][i][j2][k]*lhs[n+4][i][j][k]; // uses j2, n, updated j, n
          rhs[m][i][j2][k] = rhs[m][i][j2][k] -
            lhs[n+0][i][j2][k]*rhs[m][i][j][k]; // uses j2, n, updated j, m
        }
      }
    }
    // Forward elimination boundary case for j = grid_points[1]-2 and j1 = grid_points[1]-1 (m = 3, 4)
    j  = grid_points[1]-2;
    j1 = grid_points[1]-1;
    // These nested loops iterate over the i-k plane. Calculations for different (i, k) are independent for a fixed m.
    // This is also a target for OpenMP parallelization.
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
        // Calculations for m = 3, 4 are inside the i-k loops
        fac1               = 1./lhs[n+2][i][j][k]; // uses j, n (depends on m)
        lhs[n+3][i][j][k]   = fac1*lhs[n+3][i][j][k]; // updates j, n
        lhs[n+4][i][j][k]   = fac1*lhs[n+4][i][j][k]; // updates j, n
        rhs[m][i][j][k]     = fac1*rhs[m][i][j][k]; // updates j, m
        // Updates at j+1 (j1) use previously updated values at j
        lhs[n+2][i][j1][k] = lhs[n+2][i][j1][k] -
          lhs[n+1][i][j1][k]*lhs[n+3][i][j][k]; // uses j1, n, updated j, n
        lhs[n+3][i][j1][k] = lhs[n+3][i][j1][k] -
          lhs[n+1][i][j1][k]*lhs[n+4][i][j][k]; // uses j1, n, updated j, n
        rhs[m][i][j1][k]   = rhs[m][i][j1][k] -
          lhs[n+1][i][j1][k]*rhs[m][i][j][k]; // uses j1, n, updated j, m
        // Special handling for the last plane (j1) after update
        fac2               = 1./lhs[n+2][i][j1][k]; // uses j1, n
        rhs[m][i][j1][k]   = fac2*rhs[m][i][j1][k]; // updates j1, m
      }
    }
  }

  //-------------------------------------------------------------------------
  // Phase 2: Perform backward substitution step along the y-direction (j).
  //          This process updates rhs arrays to find the solution.
  //          The outer j-loop is sequential (decreasing) due to data dependencies.
  //          The inner i and k loops operate on independent grid points
  //          within a given set of j-planes and are parallelizable.
  //-------------------------------------------------------------------------

  // Backward substitution boundary case for j = grid_points[1]-2 and j1 = grid_points[1]-1
  j  = grid_points[1]-2;
  j1 = grid_points[1]-1;
  n = 0; // n is fixed at 0 for m = 0, 1, 2
  // Sequential loop over m
  for (m = 0; m < 3; m++) {
    // These nested loops iterate over the i-k plane. Calculations for different (i, k) are independent for a fixed m.
    // This is a target for OpenMP parallelization.
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
        // Updates at j use values from j1 (which was solved in the forward pass boundary)
	rhs[m][i][j][k] = rhs[m][i][j][k] -
	  lhs[n+3][i][j][k]*rhs[m][i][j1][k]; // uses j, j1, m, n
      }
    }
  }

  // Backward substitution boundary case for j = grid_points[1]-2 and j1 = grid_points[1]-1 (m = 3, 4)
  // Sequential loop over m
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5; // n depends on m
    // These nested loops iterate over the i-k plane. Calculations for different (i, k) are independent for a fixed m.
    // This is a target for OpenMP parallelization.
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
        // Updates at j use values from j1 (which was solved in the forward pass boundary)
	rhs[m][i][j][k] = rhs[m][i][j][k] -
	  lhs[n+3][i][j][k]*rhs[m][i][j1][k]; // uses j, j1, m, n
      }
    }
  }


  // Backward substitution main loop for j from grid_points[1]-3 down to 0 (m = 0, 1, 2)
  n = 0; // n is fixed at 0 for m = 0, 1, 2
  // Sequential loop over m
  for (m = 0; m < 3; m++) {
    // Sequential loop over j (decreasing) due to dependencies (j depends on j+1, j+2)
    for (j = grid_points[1]-3; j >= 0; j--) {
      j1 = j  + 1;
      j2 = j  + 2;
      // These nested loops iterate over the i-k plane.
      // Calculations for different (i, k) pairs are independent for a fixed m and j.
      // This is a prime target for OpenMP parallelization (e.g., #pragma omp for collapse(2))
      for (i = 1; i <= grid_points[0]-2; i++) {
        for (k = 1; k <= grid_points[2]-2; k++) {
          // Updates at j use values from j+1 (j1) and j+2 (j2) which were solved in previous j steps
          rhs[m][i][j][k] = rhs[m][i][j][k] -
            lhs[n+3][i][j][k]*rhs[m][i][j1][k] -
            lhs[n+4][i][j][k]*rhs[m][i][j2][k]; // uses j, j1, j2, m, n
        }
      }
    }
  }

  // Backward substitution main loop for j from grid_points[1]-3 down to 0 (m = 3, 4)
  // Sequential loop over m
  for (m = 3; m < 5; m++) {
    n = (m-3+1)*5; // n depends on m
    // Sequential loop over j (decreasing) due to dependencies (j depends on j+1, j+2)
    for (j = grid_points[1]-3; j >= 0; j--) {
      j1 = j  + 1;
      j2 = j1 + 1; // Same as j+2
      // These nested loops iterate over the i-k plane.
      // Calculations for different (i, k) pairs are independent for a fixed m and j.
      // This is a prime target for OpenMP parallelization (e.g., #pragma omp for collapse(2))
      for (i = 1; i <= grid_points[0]-2; i++) {
        for (k = 1; k <= grid_points[2]-2; k++) {
          // Updates at j use values from j+1 (j1) and j+2 (j2) which were solved in previous j steps
          rhs[m][i][j][k] = rhs[m][i][j][k] -
            lhs[n+3][i][j][k]*rhs[m][i][j1][k] -
            lhs[n+4][i][j][k]*rhs[m][i][j2][k]; // uses j, j1, j2, m, n
        }
      }
    }
  }

  //-------------------------------------------------------------------------
  // Phase 3: Post-processing (implicitly done by pinvr)
  //-------------------------------------------------------------------------
  pinvr();
}