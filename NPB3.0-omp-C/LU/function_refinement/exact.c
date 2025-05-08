static void exact( int i, int j, int k, double u000ijk[5] ) {
  double xi, eta, zeta;
  // Assumes nx0, ny0, nz are defined elsewhere
  xi  = ((double)i) / (nx0 - 1);
  eta  = ((double)j) / (ny0 - 1);
  zeta = ((double)k) / (nz - 1);

  // Pre-calculate polynomial terms. This calculation is common
  // for all components of u000ijk for a given i, j, k.
  // This decouples the polynomial calculation from the final sum.
  double poly_terms[13];
  poly_terms[0] = 1.0;
  poly_terms[1] = xi;
  poly_terms[2] = eta;
  poly_terms[3] = zeta;
  poly_terms[4] = xi * xi;
  poly_terms[5] = eta * eta;
  poly_terms[6] = zeta * zeta;
  poly_terms[7] = xi * xi * xi;
  poly_terms[8] = eta * eta * eta;
  poly_terms[9] = zeta * zeta * zeta;
  poly_terms[10] = xi * xi * xi * xi;
  poly_terms[11] = eta * eta * eta * eta;
  poly_terms[12] = zeta * zeta * zeta * zeta;

  // Calculate each component of u000ijk.
  // Each iteration m is independent, allowing potential parallelization
  // of this loop over m. Load is balanced per iteration.
  // Assumes ce is a global or accessible 2D array (5x13)
  for (int m = 0; m < 5; m++) {
    u000ijk[m] = 0.0; // Initialize reduction sum
    // The inner loop sums terms, a reduction operation.
    // Parallelizing this small inner loop directly is less efficient
    // than parallelizing the outer loop over m or the caller loops.
    for (int l = 0; l < 13; l++) {
        u000ijk[m] += ce[m][l] * poly_terms[l];
    }
  }
}