static void initialize(void) {
{
  int i, j, k, m, ix, iy, iz;
  double xi, eta, zeta, Pface[2][3][5], Pxi, Peta, Pzeta, temp[5];
  
  // Phase 1: Initialize all points to 1.0
  // This loop is highly parallelizable across i, j, and k dimensions.
  // Each iteration writes to a unique element u[i][j][k][m].
  for (i = 0; i < IMAX; i++) {
    for (j = 0; j < IMAX; j++) {
      for (k = 0; k < IMAX; k++) {
	for (m = 0; m < 5; m++) {
	  u[i][j][k][m] = 1.0;
	}
      }
    }
  }

  // Phase 2: Calculate interior points based on exact solution
  // This loop is highly parallelizable across i, j, and k dimensions.
  // Each iteration (i, j, k) calculates values for u[i][j][k][0..4]
  // using temporary variables (xi, eta, zeta, Pface, Pxi, Peta, Pzeta, temp)
  // that are local to this (i, j, k) iteration.
  // There are no loop-carried dependencies on u.
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      for (k = 0; k < grid_points[2]; k++) {
	zeta = (double)k * dnzm1;

	// Calculations using exact_solution are local to (i, j, k)
	for (ix = 0; ix < 2; ix++) {
	  exact_solution((double)ix, eta, zeta, 
                         &(Pface[ix][0][0]));
	}
	for (iy = 0; iy < 2; iy++) {
	  exact_solution(xi, (double)iy , zeta, 
                         &Pface[iy][1][0]);
	}
	for (iz = 0; iz < 2; iz++) {
	  exact_solution(xi, eta, (double)iz,   
                         &Pface[iz][2][0]);
	}

	// Combine results - operations are local to (i, j, k)
	for (m = 0; m < 5; m++) {
	  Pxi   = xi   * Pface[1][0][m] + 
	    (1.0-xi)   * Pface[0][0][m];
	  Peta  = eta  * Pface[1][1][m] + 
	    (1.0-eta)  * Pface[0][1][m];
	  Pzeta = zeta * Pface[1][2][m] + 
	    (1.0-zeta) * Pface[0][2][m];
	  
	  u[i][j][k][m] = Pxi + Peta + Pzeta - 
	    Pxi*Peta - Pxi*Pzeta - Peta*Pzeta + 
	    Pxi*Peta*Pzeta;
	}
      }
    }
  }

  // Phase 3: Set boundary conditions on faces.
  // These blocks are sequential in terms of faces (i=0, i=imax-1, etc.),
  // but each block is parallelizable internally over the two iterating dimensions.
  // The sequential execution of these blocks handles potential overlaps
  // at edges and corners according to the original code's order.
  
  // Boundary i=0 face: Parallelizable over j, k.
  i = 0;
  xi = 0.0;
  for (j = 0; j < grid_points[1]; j++) {
    eta = (double)j * dnym1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp); // temp is local to (j,k) iteration
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m]; // Writes to u[0][j][k][m], distinct for each (j,k,m)
      }
    }
  }

  // Boundary i=grid_points[0]-1 face: Parallelizable over j, k.
  i = grid_points[0]-1;
  xi = 1.0;
  for (j = 0; j < grid_points[1]; j++) {
    eta = (double)j * dnym1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp); // temp is local to (j,k) iteration
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m]; // Writes to u[imax-1][j][k][m], distinct for each (j,k,m)
      }
    }
  }

  // Boundary j=0 face: Parallelizable over i, k.
  j = 0;
  eta = 0.0;
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp); // temp is local to (i,k) iteration
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m]; // Writes to u[i][0][k][m], distinct for each (i,k,m)
      }
    }
  }

  // Boundary j=grid_points[1]-1 face: Parallelizable over i, k.
  j = grid_points[1]-1;
  eta = 1.0;
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp); // temp is local to (i,k) iteration
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m]; // Writes to u[i][jmax-1][k][m], distinct for each (i,k,m)
      }
    }
  }

  // Boundary k=0 face: Parallelizable over i, j.
  k = 0;
  zeta = 0.0;
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i *dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      exact_solution(xi, eta, zeta, temp); // temp is local to (i,j) iteration
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m]; // Writes to u[i][j][0][m], distinct for each (i,j,m)
      }
    }
  }

  // Boundary k=grid_points[2]-1 face: Parallelizable over i, j.
  k = grid_points[2]-1;
  zeta = 1.0;
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      exact_solution(xi, eta, zeta, temp); // temp is local to (i,j) iteration
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m]; // Writes to u[i][j][kmax-1][m], distinct for each (i,j,m)
      }
    }
  }
}
}