static void rhs_norm(double rms[5]) {
  int i, j, k, d, m;
  double add;
  
  // Initialize rms array
  #pragma omp parallel for
  for (m = 0; m < 5; m++) {
    rms[m] = 0.0;
  }
  
  // Main loop to compute RMS
  #pragma omp parallel for private(j, k, m, add) collapse(3)
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
        for (m = 0; m < 5; m++) {
          add = rhs[i][j][k][m];
          #pragma omp atomic
          rms[m] = rms[m] + add*add;
        }
      }
    }
  }

  // Normalize and take the square root of rms
  #pragma omp parallel for
  for (m = 0; m < 5; m++) {
    for (d = 0; d <= 2; d++) {
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    }
    rms[m] = sqrt(rms[m]);
  }
}