static void initialize(void) {
  int i, j, k, m, ix, iy, iz;
  double xi, eta, zeta, Pface[2][3][5], Pxi, Peta, Pzeta, temp[5];
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
  for (i = 0; i <= grid_points[0]-1; i++) {
    xi = (double)i * dnxm1;
    for (j = 0; j <= grid_points[1]-1; j++) {
      eta = (double)j * dnym1;
      for (k = 0; k <= grid_points[2]-1; k++) {
	zeta = (double)k * dnzm1;
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
	for (m = 0; m < 5; m++) {
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
  xi = 0.0;
  i  = 0;
  for (j = 0; j < grid_points[1]; j++) {
    eta = (double)j * dnym1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[m][i][j][k] = temp[m];
      }
    }
  }
  xi = 1.0;
  i  = grid_points[0]-1;
  for (j = 0; j < grid_points[1]; j++) {
    eta = (double)j * dnym1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[m][i][j][k] = temp[m];
      }
    }
  }
  eta = 0.0;
  j   = 0;
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[m][i][j][k] = temp[m];
      }
    }
  }
  eta = 1.0;
  j   = grid_points[1]-1;
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[m][i][j][k] = temp[m];
      }
    }
  }
  zeta = 0.0;
  k    = 0;
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i *dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[m][i][j][k] = temp[m];
      }
    }
  }
  zeta = 1.0;
  k    = grid_points[2]-1;
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[m][i][j][k] = temp[m];
      }
    }
  }
}