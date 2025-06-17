static void initialize(void) {
  int i, j, k, m, ix, iy, iz;
  double xi, eta, zeta, Pface[2][3][5], Pxi, Peta, Pzeta, temp[5];

  #pragma omp parallel for private(i, j, k, m) shared(u)
  for (i = 0; i < IMAX; i++) {
    for (j = 0; j < IMAX; j++) {
      for (k = 0; k < IMAX; k++) {
        for (m = 0; m < 5; m++) {
          u[i][j][k][m] = 1.0;
        }
      }
    }
  }

  #pragma omp parallel for private(i, j, k, xi, eta, zeta, ix, iy, iz, Pface, Pxi, Peta, Pzeta, m) shared(u)
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      for (k = 0; k < grid_points[2]; k++) {
        zeta = (double)k * dnzm1;
        for (ix = 0; ix < 2; ix++) {
          exact_solution((double)ix, eta, zeta, &(Pface[ix][0][0]));
        }
        for (iy = 0; iy < 2; iy++) {
          exact_solution(xi, (double)iy , zeta, &Pface[iy][1][0]);
        }
        for (iz = 0; iz < 2; iz++) {
          exact_solution(xi, eta, (double)iz, &Pface[iz][2][0]);
        }
        for (m = 0; m < 5; m++) {
          Pxi   = xi   * Pface[1][0][m] + (1.0-xi)   * Pface[0][0][m];
          Peta  = eta  * Pface[1][1][m] + (1.0-eta)  * Pface[0][1][m];
          Pzeta = zeta * Pface[1][2][m] + (1.0-zeta) * Pface[0][2][m];
          u[i][j][k][m] = Pxi + Peta + Pzeta - Pxi*Peta - Pxi*Pzeta - Peta*Pzeta + Pxi*Peta*Pzeta;
        }
      }
    }
  }

  i = 0;
  xi = 0.0;
  #pragma omp parallel for private(j, k, eta, zeta, temp, m) shared(u)
  for (j = 0; j < grid_points[1]; j++) {
    eta = (double)j * dnym1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[i][j][k][m] = temp[m];
      }
    }
  }

  i = grid_points[0]-1;
  xi = 1.0;
  #pragma omp parallel for private(j, k, eta, zeta, temp, m) shared(u)
  for (j = 0; j < grid_points[1]; j++) {
    eta = (double)j * dnym1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[i][j][k][m] = temp[m];
      }
    }
  }

  j = 0;
  eta = 0.0;
  #pragma omp parallel for private(i, k, xi, zeta, temp, m) shared(u)
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[i][j][k][m] = temp[m];
      }
    }
  }

  j = grid_points[1]-1;
  eta = 1.0;
  #pragma omp parallel for private(i, k, xi, zeta, temp, m) shared(u)
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[i][j][k][m] = temp[m];
      }
    }
  }

  k = 0;
  zeta = 0.0;
  #pragma omp parallel for private(i, j, xi, eta, temp, m) shared(u)
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[i][j][k][m] = temp[m];
      }
    }
  }

  k = grid_points[2]-1;
  zeta = 1.0;
  #pragma omp parallel for private(i, j, xi, eta, temp, m) shared(u)
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[i][j][k][m] = temp[m];
      }
    }
  }
}