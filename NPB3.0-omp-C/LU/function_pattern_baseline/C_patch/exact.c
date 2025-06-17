static void exact( int i, int j, int k, double u000ijk[5] ) {
  int m;
  double xi, eta, zeta;
  xi  = ((double)i) / (nx0 - 1);
  eta  = ((double)j) / (ny0 - 1);
  zeta = ((double)k) / (nz - 1);
#pragma omp parallel for
  for (m = 0; m < 5; m++) {
    u000ijk[m] =  ce[m][0]
      + ce[m][1] * xi
      + ce[m][2] * eta
      + ce[m][3] * zeta
      + ce[m][4] * xi * xi
      + ce[m][5] * eta * eta
      + ce[m][6] * zeta * zeta
      + ce[m][7] * xi * xi * xi
      + ce[m][8] * eta * eta * eta
      + ce[m][9] * zeta * zeta * zeta
      + ce[m][10] * xi * xi * xi * xi
      + ce[m][11] * eta * eta * eta * eta
      + ce[m][12] * zeta * zeta * zeta * zeta;
  }
}