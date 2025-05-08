// Assuming necessary headers like stdio.h and definitions for
// nx, ny, nz, nx0, ny0, u, and the exact function are available.
// For example:
// #include <stdio.h>
// #define nx ...
// #define ny ...
// #define nz ...
// #define nx0 ...
// #define ny0 ...
// double u[nx][ny][nz][5]; // Example declaration
// void exact(int x, int y, int z, double* result); // Example declaration

static void setiv_refined(void) {
  int i, j, k, m;
  int iglob, jglob; // Global indices corresponding to local i, j
  double xi, eta, zeta; // Normalized coordinates within the global grid [0, 1]
  double pxi, peta, pzeta; // Interpolated boundary contributions

  // Note: This implementation uses Variable Length Arrays (VLA), which requires C99 or later.
  // For C89 compatibility, these arrays would need to be dynamically allocated using malloc.

  // Phase 1: Pre-calculate boundary values using the exact function.
  // This phase decouples the potentially expensive exact() calls from the
  // main computational loop. The results are stored in temporary arrays.
  // Each loop nest in this phase computes values for a specific boundary plane/line.
  // These loops are independent of each other and can be parallelized.

  // Temporary storage for boundary values. Each element stores a 5-component vector.
  // Sizes are based on the index ranges required for lookups in Phase 2.
  // The k index for x/y boundaries and i/j index for z boundaries match the dimensions
  // needed by exact(ix, iy, iz).
  double boundary_x0[ny][nz][5];
  double boundary_xn1[ny][nz][5];
  double boundary_y0[nx][nz][5];
  double boundary_yn1[nx][nz][5];
  double boundary_z0[nx][ny][5];
  double boundary_zn1[nx][ny][5];

  // Compute boundary values on the global faces ix=0 and ix=nx0-1
  // These values depend on jglob and k. This loop nest is parallelizable.
  for (j = 0; j < ny; j++) {
    jglob = j; // Assuming local j maps directly to global j for exact calls
    for (k = 0; k < nz; k++) { // exact function might need k from 0 to nz-1
      exact (0, jglob, k, boundary_x0[j][k]);
      exact (nx0-1, jglob, k, boundary_xn1[j][k]);
    }
  }

  // Compute boundary values on the global faces iy=0 and iy=ny0-1
  // These values depend on iglob and k. This loop nest is parallelizable.
  for (i = 0; i < nx; i++) {
    iglob = i; // Assuming local i maps directly to global i for exact calls
    for (k = 0; k < nz; k++) { // exact function might need k from 0 to nz-1
      exact (iglob, 0, k, boundary_y0[i][k]);
      exact (iglob, ny0-1, k, boundary_yn1[i][k]);
    }
  }

  // Compute boundary values on the global faces iz=0 and iz=nz-1
  // These values depend on iglob and jglob. This loop nest is parallelizable.
  for (i = 0; i < nx; i++) {
    iglob = i; // Assuming local i maps directly to global i for exact calls
    for (j = 0; j < ny; j++) {
      jglob = j; // Assuming local j maps directly to global j for exact calls
      exact (iglob, jglob, 0, boundary_z0[i][j]);
      exact (iglob, jglob, nz-1, boundary_zn1[i][j]);
    }
  }

  // Phase 2: Compute solution points for the interior of the grid block.
  // This loop nest calculates u[i][j][k][m] for each point in the block,
  // using the pre-calculated boundary values.
  // The computation for each (i, j, k) point is independent of others within this loop nest.
  // This main loop nest is the primary target for parallelization (over i, j, or k).
  // The workload within the inner 'if' condition might be uneven due to skipping
  // local points that are on global boundaries (iglob/jglob == 0 or nx0-1/ny0-1).
  // A dynamic or guided OpenMP schedule might help with load balancing if these boundary checks are significant.

  for (j = 0; j < ny; j++) {
    jglob = j; // Assuming local j maps directly to global j for boundary checks and eta calculation
    // Skip computation for points on the global boundary planes in j
    if (jglob != 0 && jglob != ny0-1) {
      // eta is constant for a given jglob (and thus j)
      eta = ( (double) (jglob) ) / (ny0-1);

      for (k = 1; k < nz - 1; k++) { // Original k loop range (1 to nz-2), excluding global k boundaries
        // zeta is constant for a given k
        zeta = ((double)k) / (nz-1);

        // Get pointers to relevant boundary data for this (j, k) slice.
        // These values were pre-calculated in Phase 1. They are read-only here.
        double* ue_1jk = boundary_x0[j][k]; // Boundary data at global ix=0, iy=jglob, iz=k
        double* ue_nx0jk = boundary_xn1[j][k]; // Boundary data at global ix=nx0-1, iy=jglob, iz=k

        for (i = 0; i < nx; i++) {
          iglob = i; // Assuming local i maps directly to global i for boundary checks and xi calculation
          // Skip computation for points on the global boundary lines in i
          if(iglob != 0 && iglob != nx0-1) {
            // xi is constant for a given iglob (and thus i)
            xi = ( (double) (iglob) ) / (nx0-1);

            // Get pointers to relevant boundary data for this (i, j, k) point.
            // These values were pre-calculated in Phase 1. They are read-only here.
            double* ue_i1k = boundary_y0[i][k]; // Boundary data at global ix=iglob, iy=0, iz=k
            double* ue_iny0k = boundary_yn1[i][k]; // Boundary data at global ix=iglob, iy=ny0-1, iz=k
            double* ue_ij1 = boundary_z0[i][j]; // Boundary data at global ix=iglob, iy=jglob, iz=0
            double* ue_ijnz = boundary_zn1[i][j]; // Boundary data at global ix=iglob, iy=jglob, iz=nz-1

            // Compute the interpolated contributions from each boundary direction.
            // pxi, peta, pzeta are local variables for the current (i, j, k) point.
            pxi =   ( 1.0 - xi ) * ue_1jk[m] + xi   * ue_nx0jk[m];
            peta =  ( 1.0 - eta ) * ue_i1k[m] + eta   * ue_iny0k[m];
            pzeta = ( 1.0 - zeta ) * ue_ij1[m] + zeta   * ue_ijnz[m];

            // Compute the final solution components at u[i][j][k][m].
            // The loop over m computes the 5 components for the current (i, j, k) point.
            // This inner loop is very small (m=5) and usually not parallelized directly,
            // but is covered by parallelizing the outer i, j, or k loops.
            for (m = 0; m < 5; m++) {
              u[i][j][k][m] = pxi + peta + pzeta
                - pxi * peta - peta * pzeta - pzeta * pxi
                + pxi * peta * pzeta;
            } // End m loop (computes 5 components)
          } // End if iglob != 0 && iglob != nx0-1 (skip points on global i boundary)
        } // End i loop
      } // End k loop (from 1 to nz-2, skipping global k boundaries)
    } // End if jglob != 0 && jglob != ny0-1 (skip points on global j boundary)
  } // End j loop

  // Note: If using dynamic allocation for boundary arrays (malloc),
  // remember to free the allocated memory here.
  // Example for malloc:
  // double (*boundary_x0)[nz][5] = malloc(ny * sizeof(*boundary_x0));
  // ... use boundary_x0 ...
  // free(boundary_x0);
}