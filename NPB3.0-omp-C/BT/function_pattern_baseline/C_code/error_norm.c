static void error_norm(double rms[5]) {
  int i, j, k, m, d;
  double xi, eta, zeta, u_exact[5], add;

  #pragma omp parallel for private(i, j, k, xi, eta, zeta, u_exact, add) shared(rms)
  for (m = 0; m < 5; m++) {
    rms[m] = 0.0;
  }

  #pragma omp parallel for private(i, j, k, xi, eta, zeta, u_exact, add) shared(rms) collapse(3)
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      for (k = 0; k < grid_points[2]; k++) {
        zeta = (double)k * dnzm1;
        exact_solution(xi, eta, zeta, u_exact);
        for (m = 0; m < 5; m++) {
          add = u[i][j][k][m] - u_exact[m];
          #pragma omp atomic
          rms[m] = rms[m] + add*add;
        }
      }
    }
  }

  for (m = 0; m < 5; m++) {
    for (d = 0; d <= 2; d++) {
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    }
    rms[m] = sqrt(rms[m]);
  }
}