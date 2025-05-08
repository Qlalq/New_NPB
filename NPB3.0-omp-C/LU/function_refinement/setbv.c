static void setbv(void) {
  int i, j, k;
  // iglob and jglob seem to represent global coordinates in the original code,
  // calculated based on local indices i, j.
  // Explicitly calculating them per iteration makes their dependence on
  // loop variables clear, which is beneficial for parallelization.
  int iglob, jglob;

  // 1. Set boundary condition on Z faces (k=0 and k=nz-1)
  // Iterate over the i-j plane. Each (i, j) pair determines two boundary points
  // (one on k=0, one on k=nz-1).
  // Operations for different (i, j) pairs are independent.
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      // Calculate global i, j coordinates for this point.
      // Assuming local origin (0,0) corresponds to global origin (0,0)
      // based on the original code's iglob=i, jglob=j assignments.
      iglob = i; // Global i-coordinate for this point
      jglob = j; // Global j-coordinate for this point

      // Set boundary at k=0 face
      // The third argument '0' is the global k-coordinate for this face.
      exact( iglob, jglob, 0, &u[i][j][0][0] );

      // Set boundary at k=nz-1 face
      // The third argument 'nz-1' is the global k-coordinate for this face.
      exact( iglob, jglob, nz-1, &u[i][j][nz-1][0] );
    }
  }

  // 2. Set boundary condition on Y faces (j=0)
  // Iterate over the i-k plane. Each (i, k) pair determines one boundary point.
  // Operations for different (i, k) pairs are independent.
  for (i = 0; i < nx; i++) {
    for (k = 0; k < nz; k++) {
      // Calculate global i coordinate for this point.
      iglob = i; // Global i-coordinate for this point

      // Set boundary at j=0 face
      // The second argument '0' is the global j-coordinate for this face.
      // The third argument 'k' is the global k-coordinate for this point.
      exact( iglob, 0, k, &u[i][0][k][0] );
    }
  }

  // 3. Set boundary condition on Y faces (j=ny-1)
  // Iterate over the i-k plane. Each (i, k) pair determines one boundary point.
  // Operations for different (i, k) pairs are independent.
  for (i = 0; i < nx; i++) {
    for (k = 0; k < nz; k++) {
      // Calculate global i coordinate for this point.
      iglob = i; // Global i-coordinate for this point

      // Set boundary at j=ny-1 face
      // The second argument 'ny0-1' is the global j-coordinate for this face
      // based on the original code's argument.
      // The third argument 'k' is the global k-coordinate for this point.
      exact( iglob, ny0-1,  k, &u[i][ny-1][k][0] );
    }
  }

  // 4. Set boundary condition on X faces (i=0)
  // Iterate over the j-k plane. Each (j, k) pair determines one boundary point.
  // Operations for different (j, k) pairs are independent.
  for (j = 0; j < ny; j++) {
    for (k = 0; k < nz; k++) {
      // Calculate global j coordinate for this point.
      jglob = j; // Global j-coordinate for this point

      // Set boundary at i=0 face
      // The first argument '0' is the global i-coordinate for this face.
      // The second argument 'jglob' is the global j-coordinate for this point.
      // The third argument 'k' is the global k-coordinate for this point.
      exact( 0, jglob, k, &u[0][j][k][0] );
    }
  }

  // 5. Set boundary condition on X faces (i=nx-1)
  // Iterate over the j-k plane. Each (j, k) pair determines one boundary point.
  // Operations for different (j, k) pairs are independent.
  for (j = 0; j < ny; j++) {
    for (k = 0; k < nz; k++) {
      // Calculate global j coordinate for this point.
      jglob = j; // Global j-coordinate for this point

      // Set boundary at i=nx-1 face
      // The first argument 'nx0-1' is the global i-coordinate for this face
      // based on the original code's argument.
      // The second argument 'jglob' is the global j-coordinate for this point.
      // The third argument 'k' is the global k-coordinate for this point.
      exact( nx0-1, jglob, k, &u[nx-1][j][k][0] );
    }
  }
}