static void rhs_norm(double rms[5]) {
  int i, j, k, d, m;
  double add;
  for (m = 0; m < 5; m++) {
    rms[m] = 0.0;
  }
  for (i = 0; i <= grid_points[0]-2; i++) {
    for (j = 0; j <= grid_points[1]-2; j++) {
      for (k = 0; k <= grid_points[2]-2; k++) {
	for (m = 0; m < 5; m++) {
	  add = rhs[m][i][j][k];
	  rms[m] = rms[m] + add*add;
	}
      }
    }
  }
  for (m = 0; m < 5; m++) {
    for (d = 0; d < 3; d++) {
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    }
    rms[m] = sqrt(rms[m]);
  }
}