static void initialize(void) {
  int i, j, k, m, ix, iy, iz;

  // Loop 1: Initial fill (IMAX)
  // This loop initializes the entire u array buffer space.
  // It is embarrassingly parallel across i, j, and k dimensions.
  // Each iteration writes to a unique location in the u array.
  for (i = 0; i <= IMAX-1; i++) {
    for (j = 0; j <= IMAX-1; j++) {
      for (k = 0; k <= IMAX-1; k++) {
	u[0][i][j][k] = 1.0;
	u[1][i][j][k] = 0.0;
	u[2][i][j][k] = 0.0;
	u[3][i][j][k] = 0.0;
	u[4][i][j][k] = 1.0;
      }
    }
  }

  // Loop 2: Calculate Interior points using interpolation
  // The original loop iterated over 0 to grid_points[dim]-1, including boundaries.
  // To enable efficient parallelization and avoid data conflicts with the subsequent
  // boundary updates (which use the exact_solution), this loop is restricted
  // to the interior points: 1 to grid_points[dim]-2.
  // The calculation for each (i,j,k) point is independent of other (i',j',k') points.
  for (i = 1; i <= grid_points[0]-2; i++) {
    double xi = (double)i * dnxm1;
    for (j = 1; j <= grid_points[1]-2; j++) {
      double eta = (double)j * dnym1;
      for (k = 1; k <= grid_points[2]-2; k++) {
	double zeta = (double)k * dnzm1;
        // Pface is temporary storage specific to this (i,j,k) iteration.
        double Pface[2][3][5];

	// These loops calculate values from exact_solution at face locations
        // relative to the current interior point (xi, eta, zeta).
        // They are independent for a fixed (i,j,k).
	for (ix = 0; ix < 2; ix++) {
	  exact_solution((double)ix, eta, zeta,
			 &Pface[ix][0][0]);
	}
	for (iy = 0; iy < 2; iy++) {
	  exact_solution(xi, (double)iy , zeta,
			 &Pface[iy][1][0]);
	}
	for (iz = 0; iz < 2; iz++) {
	  exact_solution(xi, eta, (double)iz,
			 &Pface[iz][2][0]);
	}

        // This loop calculates the interpolated value for each component m.
        // It is independent for each m.
	for (m = 0; m < 5; m++) {
          // Pxi, Peta, Pzeta are temporary local variables for the calculation
          double Pxi, Peta, Pzeta;
	  Pxi   = xi   * Pface[1][0][m] +
	    (1.0-xi)   * Pface[0][0][m];
	  Peta  = eta  * Pface[1][1][m] +
	    (1.0-eta)  * Pface[0][1][m];
	  Pzeta = zeta * Pface[1][2][m] +
	    (1.0-zeta) * Pface[0][2][m];

	  u[m][i][j][k] = Pxi + Peta + Pzeta -
	    Pxi*Peta - Pxi*Pzeta - Peta*Pzeta +
	    Pxi*Peta*Pzeta;
	}
      }
    }
  }

  // Boundary loops: Calculate Face points using exact_solution
  // These loops set the boundary values using the exact solution.
  // The original code executed these sequentially, meaning later faces
  // (k=max, then j=max, then i=max) overwrite values on edges and corners.
  // To maintain this behavior and avoid race conditions on overlapping boundary
  // regions (edges and corners), these blocks are kept sequential.
  // However, *each loop body* (over the two varying dimensions) is parallelizable
  // internally, as each iteration writes to a unique point on that face.
  // temp[5] is temporary storage local to each iteration of the inner two loops.

  // Face i = 0
  { // Use a scope block to group related variables if needed, or for clarity
    double xi = 0.0;
    int i_const  = 0; // The constant index for this face
    for (j = 0; j < grid_points[1]; j++) {
      double eta = (double)j * dnym1;
      for (k = 0; k < grid_points[2]; k++) {
        double zeta = (double)k * dnzm1;
        // temp is local to this (j,k) iteration
        double temp[5];
        exact_solution(xi, eta, zeta, temp);
        for (m = 0; m < 5; m++) {
          u[m][i_const][j][k] = temp[m];
        }
      }
    }
  }

  // Face i = grid_points[0]-1
  {
    double xi = 1.0;
    int i_const  = grid_points[0]-1;
    for (j = 0; j < grid_points[1]; j++) {
      double eta = (double)j * dnym1;
      for (k = 0; k < grid_points[2]; k++) {
        double zeta = (double)k * dnzm1;
        double temp[5];
        exact_solution(xi, eta, zeta, temp);
        for (m = 0; m < 5; m++) {
          u[m][i_const][j][k] = temp[m];
        }
      }
    }
  }

  // Face j = 0
  {
    double eta = 0.0;
    int j_const   = 0;
    for (i = 0; i < grid_points[0]; i++) {
      double xi = (double)i * dnxm1;
      for (k = 0; k < grid_points[2]; k++) {
        double zeta = (double)k * dnzm1;
        double temp[5];
        exact_solution(xi, eta, zeta, temp);
        for (m = 0; m < 5; m++) {
          u[m][i][j_const][k] = temp[m];
        }
      }
    }
  }

  // Face j = grid_points[1]-1
  {
    double eta = 1.0;
    int j_const   = grid_points[1]-1;
    for (i = 0; i < grid_points[0]; i++) {
      double xi = (double)i * dnxm1;
      for (k = 0; k < grid_points[2]; k++) {
        double zeta = (double)k * dnzm1;
        double temp[5];
        exact_solution(xi, eta, zeta, temp);
        for (m = 0; m < 5; m++) {
          u[m][i][j_const][k] = temp[m];
        }
      }
    }
  }

  // Face k = 0
  {
    double zeta = 0.0;
    int k_const    = 0;
    for (i = 0; i < grid_points[0]; i++) {
      double xi = (double)i *dnxm1;
      for (j = 0; j < grid_points[1]; j++) {
        double eta = (double)j * dnym1;
        double temp[5];
        exact_solution(xi, eta, zeta, temp);
        for (m = 0; m < 5; m++) {
          u[m][i][j][k_const] = temp[m];
        }
      }
    }
  }

  // Face k = grid_points[2]-1
  {
    double zeta = 1.0;
    int k_const    = grid_points[2]-1;
    for (i = 0; i < grid_points[0]; i++) {
      double xi = (double)i * dnxm1;
      for (j = 0; j < grid_points[1]; j++) {
        double eta = (double)j * dnym1;
        double temp[5];
        exact_solution(xi, eta, zeta, temp);
        for (m = 0; m < 5; m++) {
          u[m][i][j][k_const] = temp[m];
        }
      }
    }
  }
}