// Note: This code assumes the necessary declarations for arrays
// rsd, frct, flux, and constants like nx, ny, nz, nx0, ny0,
// ist, iend, jst, jend, dssp, ce, tx2, tx3, ty2, ty3, tz2, tz3,
// C1, C2, C3, C4, C5, dx1..dx5, dy1..dy5, dz1..dz5 are available
// in the scope where this function is defined (e.g., globals or includes).

static void erhs(void) {
  // Local variables used within the function. These variables are typically
  // private to each thread when loops are parallelized, avoiding dependencies.
  int i, j, k, m;
  int iglob, jglob;
  int L1, L2; // Generic loop bounds, redefined in different sections
  int ist1, iend1; // Loop bounds for boundary regions in i-direction
  int jst1, jend1; // Loop bounds for boundary regions in j-direction
  double  dsspm;
  double  xi, eta, zeta; // Coordinates used in rsd initialization
  double  q; // Intermediate calculation for pressure term
  // Temporary variables used in flux calculations within directional passes
  double  u21, u31, u41;
  double  tmp;
  double  u21i, u31i, u41i, u51i;
  double  u21im1, u31im1, u41im1, u51im1;
  double  u21j, u31j, u41j, u51j;
  double  u21jm1, u31jm1, u41jm1, u51jm1;
  double  u21k, u31k, u41k, u51k;
  double  u21km1, u31km1, u41km1, u51km1;

  dsspm = dssp; // Assigning global dssp to local dsspm

  //---------------------------------------------------------------------
  // initialization of frct (right hand side)
  // This section initializes the frct array to zero. Each element
  // frct[i][j][k][m] is written independently.
  // This loop nest is highly parallelizable over any of the outer indices (i, j, k).
  // Parallelizing the outermost loop will distribute large chunks of work.
  // No data dependencies or synchronization needed within this section.
  //---------------------------------------------------------------------
  for (i = 0; i < nx; i++) { // Potential OpenMP parallel loop
    for (j = 0; j < ny; j++) { // Potential OpenMP parallel loop (if i not parallelized)
      for (k = 0; k < nz; k++) { // Potential OpenMP parallel loop (if i, j not parallelized)
	for (m = 0; m < 5; m++) {
	  frct[i][j][k][m] = 0.0;
	}
      }
    }
  }

  //---------------------------------------------------------------------
  // computation of the right hand side based on the residual (rsd)
  // This section computes the rsd array based on grid coordinates and constants.
  // Each element rsd[i][j][k][m] is computed independently.
  // This loop nest is highly parallelizable over any of the outer indices (i, j, k).
  // No data dependencies or synchronization needed within this section.
  //---------------------------------------------------------------------
  for (i = 0; i < nx; i++) { // Potential OpenMP parallel loop
    iglob = i;
    xi = ( (double)(iglob) ) / ( nx0 - 1 );
    for (j = 0; j < ny; j++) { // Potential OpenMP parallel loop (if i not parallelized)
      jglob = j;
      eta = ( (double)(jglob) ) / ( ny0 - 1 );
      for (k = 0; k < nz; k++) { // Potential OpenMP parallel loop (if i, j not parallelized)
	zeta = ( (double)(k) ) / ( nz - 1 );
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] =  ce[m][0]
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
    }
  }

  //---------------------------------------------------------------------
  // x-direction flux and frct calculations
  // This section processes the grid primarily along the x-direction.
  // The computations involve stencils and updates to frct and flux.
  // The calculations for different (j, k) planes are independent of each other.
  // Parallelization is most effective over the outer loops (j or k).
  // The inner loops over i within this section have spatial dependencies
  // (e.g., i-1, i+1) that generally require sequential execution for a fixed (j, k) column.
  // frct is accumulated, but writes to frct[i][j][k] are distinct for
  // different parallel iterations of the outer loop (j or k).
  // The 'flux' array is temporary within this directional pass.
  //---------------------------------------------------------------------

  // Calculate initial flux in x-direction based on rsd at the same point.
  // This calculation for flux[i][j][k] is independent for each (i,j,k).
  // Parallelize over i, j, or k.
  L1 = 0; L2 = nx-1; // Range [0, nx-1] for i
  for (i = L1; i <= L2; i++) { // Potential OpenMP parallel loop
    for (j = jst; j <= jend; j++) { // Potential OpenMP parallel loop (if i not parallelized)
      for (k = 1; k < nz - 1; k++) { // Potential OpenMP parallel loop (if i, j not parallelized)
	// Calculation for flux[i][j][k] based on rsd[i][j][k]
	flux[i][j][k][0] = rsd[i][j][k][1];
	u21 = rsd[i][j][k][1] / rsd[i][j][k][0];
	q = 0.50 * (  rsd[i][j][k][1] * rsd[i][j][k][1]
		      + rsd[i][j][k][2] * rsd[i][j][k][2]
		      + rsd[i][j][k][3] * rsd[i][j][k][3] )
	  / rsd[i][j][k][0];
	flux[i][j][k][1] = rsd[i][j][k][1] * u21 + C2 *
	  ( rsd[i][j][k][4] - q );
	flux[i][j][k][2] = rsd[i][j][k][2] * u21;
	flux[i][j][k][3] = rsd[i][j][k][3] * u21;
	flux[i][j][k][4] = ( C1 * rsd[i][j][k][4] - C2 * q ) * u21;
      }
    }
  }

  // Apply x-direction differences and viscous terms to frct.
  // This section contains multiple sequential loops over i for fixed j, k.
  // Parallelize the outer loop over j or k to process independent (j, k) columns.
  for (j = jst; j <= jend; j++) { // Recommended OpenMP parallel loop for x-direction pass
    for (k = 1; k <= nz - 2; k++) { // Potential OpenMP parallel loop (if j not parallelized, or use collapse(2))
      // Loop 4.1: Update frct using centered flux differences (i+1, i-1). Reads from flux calculated previously.
      // Writes to frct[i][j][k]. Different i iterations write to different frct[i][j][k].
      for (i = ist; i <= iend; i++) {
	for (m = 0; m < 5; m++) {
	  frct[i][j][k][m] =  frct[i][j][k][m]
	    - tx2 * ( flux[i+1][j][k][m] - flux[i-1][j][k][m] );
	}
      }
      // Loop 4.2: Recalculate flux using upwind-like differences (i, i-1) on rsd components.
      // Writes to flux[i][j][k]. Different i iterations write to different flux[i][j][k].
      for (i = ist; i <= nx-1; i++) { // L2 is nx-1
	tmp = 1.0 / rsd[i][j][k][0];
	u21i = tmp * rsd[i][j][k][1]; u31i = tmp * rsd[i][j][k][2]; u41i = tmp * rsd[i][j][k][3]; u51i = tmp * rsd[i][j][k][4];
	tmp = 1.0 / rsd[i-1][j][k][0];
	u21im1 = tmp * rsd[i-1][j][k][1]; u31im1 = tmp * rsd[i-1][j][k][2]; u41im1 = tmp * rsd[i-1][j][k][3]; u51im1 = tmp * rsd[i-1][j][k][4];

	flux[i][j][k][1] = (4.0/3.0) * tx3 * ( u21i - u21im1 );
	flux[i][j][k][2] = tx3 * ( u31i - u31im1 );
	flux[i][j][k][3] = tx3 * ( u41i - u41im1 );
	flux[i][j][k][4] = 0.50 * ( 1.0 - C1*C5 ) * tx3 *
                           ( ( u21i*u21i + u31i*u31i + u41i*u41i ) - ( u21im1*u21im1 + u31im1*u31im1 + u41im1*u41im1 ) )
	                 + (1.0/6.0) * tx3 * ( u21i*u21i - u21im1*u21im1 )
	                 + C1 * C5 * tx3 * ( u51i - u51im1 );
      }
      // Loop 4.3: Update frct using recalculated flux differences (i+1, i) and centered rsd differences (i-1, i, i+1).
      // Reads from flux recalculated in Loop 4.2. Writes to frct[i][j][k].
      // Different i iterations write to different frct[i][j][k].
      for (i = ist; i <= iend; i++) {
	frct[i][j][k][0] = frct[i][j][k][0] + dx1 * tx1 * ( rsd[i-1][j][k][0] - 2.0 * rsd[i][j][k][0] + rsd[i+1][j][k][0] );
	frct[i][j][k][1] = frct[i][j][k][1] + tx3 * C3 * C4 * ( flux[i+1][j][k][1] - flux[i][j][k][1] )
	                                  + dx2 * tx1 * ( rsd[i-1][j][k][1] - 2.0 * rsd[i][j][k][1] + rsd[i+1][j][k][1] );
	frct[i][j][k][2] = frct[i][j][k][2] + tx3 * C3 * C4 * ( flux[i+1][j][k][2] - flux[i][j][k][2] )
	                                  + dx3 * tx1 * ( rsd[i-1][j][k][2] - 2.0 * rsd[i][j][k][2] + rsd[i+1][j][k][2] );
	frct[i][j][k][3] = frct[i][j][k][3] + tx3 * C3 * C4 * ( flux[i+1][j][k][3] - flux[i][j][k][3] )
	                                  + dx4 * tx1 * ( rsd[i-1][j][k][3] - 2.0 * rsd[i][j][k][3] + rsd[i+1][j][k][3] );
	frct[i][j][k][4] = frct[i][j][k][4] + tx3 * C3 * C4 * ( flux[i+1][j][k][4] - flux[i][j][k][4] )
	                                  + dx5 * tx1 * ( rsd[i-1][j][k][4] - 2.0 * rsd[i][j][k][4] + rsd[i+1][j][k][4] );
      }
      // Loops 4.4, 4.5, 4.6: Apply higher-order (5-point) stencil for viscous terms on boundaries and interior.
      // These apply only to specific i indices (1, 2, 3..nx-4, nx-3, nx-2).
      // Writes to frct[i][j][k]. Different i indices are written to.
      // These are independent for different j and k.
      for (m = 0; m < 5; m++) {
	frct[1][j][k][m] = frct[1][j][k][m] - dsspm * ( + 5.0 * rsd[1][j][k][m] - 4.0 * rsd[2][j][k][m] + rsd[3][j][k][m] );
	frct[2][j][k][m] = frct[2][j][k][m] - dsspm * ( - 4.0 * rsd[1][j][k][m] + 6.0 * rsd[2][j][k][m] - 4.0 * rsd[3][j][k][m] + rsd[4][j][k][m] );
      }
      ist1 = 3; // Loop 4.5 range
      iend1 = nx - 4;
      for (i = ist1; i <=iend1; i++) {
	for (m = 0; m < 5; m++) {
	  frct[i][j][k][m] = frct[i][j][k][m] - dsspm * ( rsd[i-2][j][k][m] - 4.0 * rsd[i-1][j][k][m] + 6.0 * rsd[i][j][k][m] - 4.0 * rsd[i+1][j][k][m] + rsd[i+2][j][k][m] );
	}
      }
      for (m = 0; m < 5; m++) {
	frct[nx-3][j][k][m] = frct[nx-3][j][k][m] - dsspm * ( rsd[nx-5][j][k][m] - 4.0 * rsd[nx-4][j][k][m] + 6.0 * rsd[nx-3][j][k][m] - 4.0 * rsd[nx-2][j][k][m]  );
	frct[nx-2][j][k][m] = frct[nx-2][j][k][m] - dsspm * ( rsd[nx-4][j][k][m] - 4.0 * rsd[nx-3][j][k][m] + 5.0 * rsd[nx-2][j][k][m] );
      }
    } // End loop k
  } // End loop j

  //---------------------------------------------------------------------
  // y-direction flux and frct calculations
  // This section processes the grid primarily along the y-direction.
  // The computations involve stencils and updates to frct and flux.
  // The calculations for different (i, k) planes are independent of each other.
  // Parallelization is most effective over the outer loops (i or k).
  // The inner loops over j within this section have spatial dependencies
  // that generally require sequential execution for a fixed (i, k) column.
  // frct is accumulated, but writes to frct[i][j][k] are distinct for
  // different parallel iterations of the outer loop (i or k).
  // The 'flux' array is temporary within this directional pass.
  //---------------------------------------------------------------------

  // Calculate initial flux in y-direction based on rsd at the same point.
  // This calculation for flux[i][j][k] is independent for each (i,j,k).
  // Parallelize over i, j, or k.
  L1 = 0; L2 = ny-1; // Range [0, ny-1] for j
  for (i = ist; i <= iend; i++) { // Potential OpenMP parallel loop
    for (j = L1; j <= L2; j++) { // Potential OpenMP parallel loop (if i not parallelized)
      for (k = 1; k <= nz - 2; k++) { // Potential OpenMP parallel loop (if i, j not parallelized)
	// Calculation for flux[i][j][k] based on rsd[i][j][k]
	flux[i][j][k][0] = rsd[i][j][k][2]; // Note: Uses rsd[...][2] component
	u31 = rsd[i][j][k][2] / rsd[i][j][k][0];
	q = 0.50 * (  rsd[i][j][k][1] * rsd[i][j][k][1]
		      + rsd[i][j][k][2] * rsd[i][j][k][2]
		      + rsd[i][j][k][3] * rsd[i][j][k][3] )
	  / rsd[i][j][k][0];
	flux[i][j][k][1] = rsd[i][j][k][1] * u31;
	flux[i][j][k][2] = rsd[i][j][k][2] * u31 + C2 *
	  ( rsd[i][j][k][4] - q );
	flux[i][j][k][3] = rsd[i][j][k][3] * u31;
	flux[i][j][k][4] = ( C1 * rsd[i][j][k][4] - C2 * q ) * u31;
      }
    }
  }

  // Apply y-direction differences and viscous terms to frct.
  // This section contains multiple sequential loops over j for fixed i, k.
  // Parallelize the outer loop over i or k to process independent (i, k) columns.
  for (i = ist; i <= iend; i++) { // Recommended OpenMP parallel loop for y-direction pass
    for (k = 1; k <= nz - 2; k++) { // Potential OpenMP parallel loop (if i not parallelized, or use collapse(2))
      // Loop 6.1: Update frct using centered flux differences (j+1, j-1). Reads from flux calculated previously.
      // Writes to frct[i][j][k]. Different j iterations write to different frct[i][j][k].
      for (j = jst; j <= jend; j++) {
	for (m = 0; m < 5; m++) {
	  frct[i][j][k][m] =  frct[i][j][k][m]
	    - ty2 * ( flux[i][j+1][k][m] - flux[i][j-1][k][m] );
	}
      }
      // Loop 6.2: Recalculate flux using upwind-like differences (j, j-1) on rsd components.
      // Writes to flux[i][j][k]. Different j iterations write to different flux[i][j][k].
      for (j = jst; j <= ny-1; j++) { // L2 is ny-1
	tmp = 1.0 / rsd[i][j][k][0];
	u21j = tmp * rsd[i][j][k][1]; u31j = tmp * rsd[i][j][k][2]; u41j = tmp * rsd[i][j][k][3]; u51j = tmp * rsd[i][j][k][4];
	tmp = 1.0 / rsd[i][j-1][k][0];
	u21jm1 = tmp * rsd[i][j-1][k][1]; u31jm1 = tmp * rsd[i][j-1][k][2]; u41jm1 = tmp * rsd[i][j-1][k][3]; u51jm1 = tmp * rsd[i][j-1][k][4];

	flux[i][j][k][1] = ty3 * ( u21j - u21jm1 );
	flux[i][j][k][2] = (4.0/3.0) * ty3 * ( u31j - u31jm1 );
	flux[i][j][k][3] = ty3 * ( u41j - u41jm1 );
	flux[i][j][k][4] = 0.50 * ( 1.0 - C1*C5 ) * ty3 *
                           ( ( u21j  *u21j + u31j  *u31j + u41j  *u41j ) - ( u21jm1*u21jm1 + u31jm1*u31jm1 + u41jm1*u41jm1 ) )
	                 + (1.0/6.0) * ty3 * ( u31j*u31j - u31jm1*u31jm1 )
	                 + C1 * C5 * ty3 * ( u51j - u51jm1 );
      }
      // Loop 6.3: Update frct using recalculated flux differences (j+1, j) and centered rsd differences (j-1, j, j+1).
      // Reads from flux recalculated in Loop 6.2. Writes to frct[i][j][k].
      // Different j iterations write to different frct[i][j][k].
      for (j = jst; j <= jend; j++) {
	frct[i][j][k][0] = frct[i][j][k][0] + dy1 * ty1 * ( rsd[i][j-1][k][0] - 2.0 * rsd[i][j][k][0] + rsd[i][j+1][k][0] );
	frct[i][j][k][1] = frct[i][j][k][1] + ty3 * C3 * C4 * ( flux[i][j+1][k][1] - flux[i][j][k][1] )
	                                  + dy2 * ty1 * ( rsd[i][j-1][k][1] - 2.0 * rsd[i][j][k][1] + rsd[i][j+1][k][1] );
	frct[i][j][k][2] = frct[i][j][k][2] + ty3 * C3 * C4 * ( flux[i][j+1][k][2] - flux[i][j][k][2] )
	                                  + dy3 * ty1 * ( rsd[i][j-1][k][2] - 2.0 * rsd[i][j][k][2] + rsd[i][j+1][k][2] );
	frct[i][j][k][3] = frct[i][j][k][3] + ty3 * C3 * C4 * ( flux[i][j+1][k][3] - flux[i][j][k][3] )
	                                  + dy4 * ty1 * ( rsd[i][j-1][k][3] - 2.0 * rsd[i][j][k][3] + rsd[i][j+1][k][3] );
	frct[i][j][k][4] = frct[i][j][k][4] + ty3 * C3 * C4 * ( flux[i][j+1][k][4] - flux[i][j][k][4] )
	                                  + dy5 * ty1 * ( rsd[i][j-1][k][4] - 2.0 * rsd[i][j][k][4] + rsd[i][j+1][k][4] );
      }
      // Loops 6.4, 6.5, 6.6: Apply higher-order (5-point) stencil for viscous terms on boundaries and interior.
      // These apply only to specific j indices (1, 2, 3..ny-4, ny-3, ny-2).
      // Writes to frct[i][j][k]. Different j indices are written to.
      // These are independent for different i and k.
      for (m = 0; m < 5; m++) {
	frct[i][1][k][m] = frct[i][1][k][m] - dsspm * ( + 5.0 * rsd[i][1][k][m] - 4.0 * rsd[i][2][k][m] + rsd[i][3][k][m] );
	frct[i][2][k][m] = frct[i][2][k][m] - dsspm * ( - 4.0 * rsd[i][1][k][m] + 6.0 * rsd[i][2][k][m] - 4.0 * rsd[i][3][k][m] + rsd[i][4][k][m] );
      }
      jst1 = 3; // Loop 6.5 range
      jend1 = ny - 4;
      for (j = jst1; j <= jend1; j++) {
	for (m = 0; m < 5; m++) {
	  frct[i][j][k][m] = frct[i][j][k][m] - dsspm * ( rsd[i][j-2][k][m] - 4.0 * rsd[i][j-1][k][m] + 6.0 * rsd[i][j][k][m] - 4.0 * rsd[i][j+1][k][m] + rsd[i][j+2][k][m] );
	}
      }
      for (m = 0; m < 5; m++) {
	frct[i][ny-3][k][m] = frct[i][ny-3][k][m] - dsspm * ( rsd[i][ny-5][k][m] - 4.0 * rsd[i][ny-4][k][m] + 6.0 * rsd[i][ny-3][k][m] - 4.0 * rsd[i][ny-2][k][m]  );
	frct[i][ny-2][k][m] = frct[i][ny-2][k][m] - dsspm * ( rsd[i][ny-4][k][m] - 4.0 * rsd[i][ny-3][k][m] + 5.0 * rsd[i][ny-2][k][m]  );
      }
    } // End loop k
  } // End loop i

  //---------------------------------------------------------------------
  // z-direction flux and frct calculations
  // This section processes the grid primarily along the z-direction.
  // The computations involve stencils and updates to frct and flux.
  // The calculations for different (i, j) planes are independent of each other.
  // Parallelization is most effective over the outer loops (i or j).
  // The inner loops over k within this section have spatial dependencies
  // that generally require sequential execution for a fixed (i, j) column.
  // frct is accumulated, but writes to frct[i][j][k] are distinct for
  // different parallel iterations of the outer loop (i or j).
  // The 'flux' array is temporary within this directional pass.
  //---------------------------------------------------------------------

  // This block contains both initial flux calculation and frct updates / flux recalculation, nested together.
  // The outer loops over i and j are independent. Parallelize i or j to process independent (i, j) columns.
  for (i = ist; i <= iend; i++) { // Recommended OpenMP parallel loop for z-direction pass
    for (j = jst; j <= jend; j++) { // Potential OpenMP parallel loop (if i not parallelized, or use collapse(2))
      // Block 7: Calculate initial flux in z-direction for this (i, j) column.
      // Writes to flux[i][j][k] independently for each k within this (i,j).
      // This k-loop can be parallelized, potentially nested within an (i,j) parallel region.
      for (k = 0; k <= nz-1; k++) { // Potential OpenMP parallel loop (nested within i,j parallel region)
	// Calculation for flux[i][j][k] based on rsd[i][j][k]
	flux[i][j][k][0] = rsd[i][j][k][3]; // Note: Uses rsd[...][3] component
	u41 = rsd[i][j][k][3] / rsd[i][j][k][0];
	q = 0.50 * (  rsd[i][j][k][1] * rsd[i][j][k][1]
		      + rsd[i][j][k][2] * rsd[i][j][k][2]
		      + rsd[i][j][k][3] * rsd[i][j][k][3] )
	  / rsd[i][j][k][0];
	flux[i][j][k][1] = rsd[i][j][k][1] * u41;
	flux[i][j][k][2] = rsd[i][j][k][2] * u41;
	flux[i][j][k][3] = rsd[i][j][k][3] * u41 + C2 *
	  ( rsd[i][j][k][4] - q );
	flux[i][j][k][4] = ( C1 * rsd[i][j][k][4] - C2 * q ) * u41;
      }
      // Block 8.1: Update frct using centered flux differences (k+1, k-1) from Block 7.
      // Reads from flux calculated just above for this (i,j). Writes to frct[i][j][k].
      // Different k iterations write to different frct[i][j][k].
      for (k = 1; k <= nz - 2; k++) {
	for (m = 0; m < 5; m++) {
	  frct[i][j][k][m] =  frct[i][j][k][m]
	    - tz2 * ( flux[i][j][k+1][m] - flux[i][j][k-1][m] );
	}
      }
      // Block 8.2: Recalculate flux using upwind-like differences (k, k-1) on rsd components for this (i,j) column.
      // Writes to flux[i][j][k]. Different k iterations write to different flux[i][j][k].
      for (k = 1; k <= nz-1; k++) {
	tmp = 1.0 / rsd[i][j][k][0];
	u21k = tmp * rsd[i][j][k][1]; u31k = tmp * rsd[i][j][k][2]; u41k = tmp * rsd[i][j][k][3]; u51k = tmp * rsd[i][j][k][4];
	tmp = 1.0 / rsd[i][j][k-1][0];
	u21km1 = tmp * rsd[i][j][k-1][1]; u31km1 = tmp * rsd[i][j][k-1][2]; u41km1 = tmp * rsd[i][j][k-1][3]; u51km1 = tmp * rsd[i][j][k-1][4];

	flux[i][j][k][1] = tz3 * ( u21k - u21km1 );
	flux[i][j][k][2] = tz3 * ( u31k - u31km1 );
	flux[i][j][k][3] = (4.0/3.0) * tz3 * ( u41k - u41km1 );
	flux[i][j][k][4] = 0.50 * ( 1.0 - C1*C5 ) * tz3 *
                           ( ( u21k  *u21k + u31k  *u31k + u41k  *u41k ) - ( u21km1*u21km1 + u31km1*u31km1 + u41km1*u41km1 ) )
	                 + (1.0/6.0) * tz3 * ( u41k*u41k - u41km1*u41km1 )
	                 + C1 * C5 * tz3 * ( u51k - u51km1 );
      }
      // Block 8.3: Update frct using recalculated flux differences (k+1, k) and centered rsd differences (k-1, k, k+1).
      // Reads from flux recalculated in Block 8.2. Writes to frct[i][j][k].
      // Different k iterations write to different frct[i][j][k].
      for (k = 1; k <= nz - 2; k++) {
	frct[i][j][k][0] = frct[i][j][k][0] + dz1 * tz1 * ( rsd[i][j][k+1][0] - 2.0 * rsd[i][j][k][0] + rsd[i][j][k-1][0] );
	frct[i][j][k][1] = frct[i][j][k][1] + tz3 * C3 * C4 * ( flux[i][j][k+1][1] - flux[i][j][k][1] )
	                                  + dz2 * tz1 * ( rsd[i][j][k+1][1] - 2.0 * rsd[i][j][k][1] + rsd[i][j][k-1][1] );
	frct[i][j][k][2] = frct[i][j][k][2] + tz3 * C3 * C4 * ( flux[i][j][k+1][2] - flux[i][j][k][2] )
	                                  + dz3 * tz1 * ( rsd[i][j][k+1][2] - 2.0 * rsd[i][j][k][2] + rsd[i][j][k-1][2] );
	frct[i][j][k][3] = frct[i][j][k][3] + tz3 * C3 * C4 * ( flux[i][j][k+1][3] - flux[i][j][k][3] )
	                                  + dz4 * tz1 * ( rsd[i][j][k+1][3] - 2.0 * rsd[i][j][k][3] + rsd[i][j][k-1][3] );
	frct[i][j][k][4] = frct[i][j][k][4] + tz3 * C3 * C4 * ( flux[i][j][k+1][4] - flux[i][j][k][4] )
	                                  + dz5 * tz1 * ( rsd[i][j][k+1][4] - 2.0 * rsd[i][j][k][4] + rsd[i][j][k-1][4] );
      }
      // Loops 8.4, 8.5, 8.6: Apply higher-order (5-point) stencil for viscous terms on boundaries and interior.
      // These apply only to specific k indices (1, 2, 3..nz-4, nz-3, nz-2).
      // Writes to frct[i][j][k]. Different k indices are written to.
      // These are independent for different i and j.
      for (m = 0; m < 5; m++) {
	frct[i][j][1][m] = frct[i][j][1][m] - dsspm * ( + 5.0 * rsd[i][j][1][m] - 4.0 * rsd[i][j][2][m] + rsd[i][j][3][m] );
	frct[i][j][2][m] = frct[i][j][2][m] - dsspm * (- 4.0 * rsd[i][j][1][m] + 6.0 * rsd[i][j][2][m] - 4.0 * rsd[i][j][3][m] + rsd[i][j][4][m] );
      }
      // Loop 8.5 range
      for (k = 3; k <= nz - 4; k++) {
	for (m = 0; m < 5; m++) {
	  frct[i][j][k][m] = frct[i][j][k][m] - dsspm * ( rsd[i][j][k-2][m] - 4.0 * rsd[i][j][k-1][m] + 6.0 * rsd[i][j][k][m] - 4.0 * rsd[i][j][k+1][m] + rsd[i][j][k+2][m] );
	}
      }
      for (m = 0; m < 5; m++) {
	frct[i][j][nz-3][m] = frct[i][j][nz-3][m] - dsspm * ( rsd[i][j][nz-5][m] - 4.0 * rsd[i][j][nz-4][m] + 6.0 * rsd[i][j][nz-3][m] - 4.0 * rsd[i][j][nz-2][m]  );
        frct[i][j][nz-2][m] = frct[i][j][nz-2][m] - dsspm * ( rsd[i][j][nz-4][m] - 4.0 * rsd[i][j][nz-3][m] + 5.0 * rsd[i][j][nz-2][m]  );
      }
    } // End loop j (outer for Z-pass)
  } // End loop i (outer for Z-pass)
}