static void rhs(void) {
  int i, j, k, m;
  int L1, L2; // Variables used for loop bounds in the original code
  int ist1, iend1; // Variables used for inner i-loop bounds
  int jst1, jend1; // Variables used for inner j-loop bounds
  // int kst1, kend1; // Corresponding variables for k-loop (not used in original input bounds)
  double  q;
  double  u21, u31, u41; // Local variables for flux calculations (convective)
  double  tmp; // Temporary variable for calculations
  double  u21i, u31i, u41i, u51i; // Local variables for x-viscous flux calculation
  double  u21j, u31j, u41j, u51j; // Local variables for y-viscous flux calculation
  double  u21k, u31k, u41k, u51k; // Local variables for z-viscous flux calculation
  double  u21im1, u31im1, u41im1, u51im1; // Local variables for x-viscous flux calculation (i-1)
  double  u21jm1, u31jm1, u41jm1, u51jm1; // Local variables for y-viscous flux calculation (j-1)
  double  u21km1, u31km1, u41km1, u51km1; // Local variables for z-viscous flux calculation (k-1)

  //---------------------------------------------------------------------
  // Phase 1: Initialize rsd array
  // rsd[i][j][k][m] = -frct[i][j][k][m]
  // This computation is fully independent for each element (i, j, k, m).
  // Highly parallelizable across i, j, and k.
  // Data Dependencies: Reads frct, Writes rsd. No loop-carried dependencies.
  // Load Balance: Good if iterations are distributed evenly across threads.
  // Synchronization: Not required within this block.
  //---------------------------------------------------------------------
  for (i = 0; i <= nx-1; i++) {
    for (j = 0; j <= ny-1; j++) {
      for (k = 0; k <= nz-1; k++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] = - frct[i][j][k][m];
	}
      }
    }
  }

  //---------------------------------------------------------------------
  // Phase 2: X-direction calculations
  // Computes fluxes in the x-direction and updates rsd based on them.
  // This phase requires u and frct (already used) and writes to flux and rsd.
  // The reuse of the 'flux' array for Y and Z phases introduces a necessary
  // sequential dependency between Phase 2, Phase 3, and Phase 4.
  // Inside this phase, many computations are parallelizable.
  //---------------------------------------------------------------------

  // Sub-phase 2.1: Calculate convective flux based on u (x-direction)
  // flux[i][j][k][m] depends only on u[i][j][k][m] at the same index.
  // Parallelizable across i, j, and k. Loop structure is i, j, k.
  // Data Dependencies: Reads u, Writes flux. No loop-carried dependencies on flux or u.
  // Load Balance: Good if iterations are distributed evenly.
  // Synchronization: Not required within this sub-phase.
  // (L1 and L2 are used as in original code, effectively 0 to nx-1 for i)
  L1 = 0; L2 = nx-1;
  for (i = L1; i <= L2; i++) {
    for (j = jst; j <= jend; j++) {
      for (k = 1; k <= nz - 2; k++) {
	flux[i][j][k][0] = u[i][j][k][1];
	u21 = u[i][j][k][1] / u[i][j][k][0];
	q = 0.50 * (  u[i][j][k][1] * u[i][j][k][1]
		      + u[i][j][k][2] * u[i][j][k][2]
		      + u[i][j][k][3] * u[i][j][k][3] )
	  / u[i][j][k][0];
	flux[i][j][k][1] = u[i][j][k][1] * u21 + C2 *
	  ( u[i][j][k][4] - q );
	flux[i][j][k][2] = u[i][j][k][2] * u21;
	flux[i][j][k][3] = u[i][j][k][3] * u21;
	flux[i][j][k][4] = ( C1 * u[i][j][k][4] - C2 * q ) * u21;
      }
    }
  }

  // Sub-phase 2.2: Apply x-direction contributions to rsd
  // This block contains several loop nests that update rsd based on x-direction fluxes and diffusion.
  // All updates within this block for a specific (i, j, k, m) read the existing rsd[i][j][k][m]
  // and add a new term. The order of these additions within the X, Y, Z phases is fixed.
  // The loops here are structured with j, k outermost, then i, m inner.
  // Parallelizable across j and k for the outer loops, and across i for the inner loops.
  // Data Dependencies: Reads rsd, u, flux. Writes rsd.
  // Updates to rsd[i][j][k][m] depend on flux[i-1/i+1][j][k][m] or flux[i/i+1][j][k][m]
  // and u[i-2..i+2][j][k][m]. No loop-carried dependencies on rsd or flux within the i loops.
  // Load Balance: Can be good depending on how loops (j, k, i) are parallelized.
  // Synchronization: Not required within this sub-phase if parallelized carefully (e.g., different threads update different (j,k) slices or different i ranges). Implicit sync needed before Phase 3 due to 'flux' reuse.
  for (j = jst; j <= jend; j++) {
    for (k = 1; k <= nz - 2; k++) {

      // rsd update based on convective flux differences in x-direction
      // Reads flux calculated in Sub-phase 2.1. Parallelizable across i.
      for (i = ist; i <= iend; i++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] =  rsd[i][j][k][m]
	    - tx2 * ( flux[i+1][j][k][m] - flux[i-1][j][k][m] );
	}
      }

      // Calculate viscous flux in x-direction (overwrites previous flux)
      // Depends on u[i] and u[i-1]. Parallelizable across i.
      // (L2 is used as in original code, effectively nx-1 for i)
      L2 = nx-1;
      for (i = ist; i <= L2; i++) {
	tmp = 1.0 / u[i][j][k][0];
	u21i = tmp * u[i][j][k][1]; u31i = tmp * u[i][j][k][2];
	u41i = tmp * u[i][j][k][3]; u51i = tmp * u[i][j][k][4];
	tmp = 1.0 / u[i-1][j][k][0];
	u21im1 = tmp * u[i-1][j][k][1]; u31im1 = tmp * u[i-1][j][k][2];
	u41im1 = tmp * u[i-1][j][k][3]; u51im1 = tmp * u[i-1][j][k][4];
	flux[i][j][k][1] = (4.0/3.0) * tx3 * (u21i-u21im1);
	flux[i][j][k][2] = tx3 * ( u31i - u31im1 );
	flux[i][j][k][3] = tx3 * ( u41i - u41im1 );
	// Assuming pow2 is defined elsewhere (e.g., #define pow2(a) ((a)*(a)))
	flux[i][j][k][4] = 0.50 * ( 1.0 - C1*C5 ) * tx3 * ( ( pow2(u21i) + pow2(u31i) + pow2(u41i) ) - ( pow2(u21im1) + pow2(u31im1) + pow2(u41im1) ) ) + (1.0/6.0) * tx3 * ( pow2(u21i) - pow2(u21im1) ) + C1 * C5 * tx3 * ( u51i - u51im1 );
      }

      // rsd update based on diffusion terms in x-direction and viscous flux
      // Uses flux calculated in the previous loop. Parallelizable across i.
      for (i = ist; i <= iend; i++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][0] = rsd[i][j][k][0] + dx1 * tx1 * ( u[i-1][j][k][0] - 2.0 * u[i][j][k][0] + u[i+1][j][k][0] );
	  rsd[i][j][k][1] = rsd[i][j][k][1] + tx3 * C3 * C4 * ( flux[i+1][j][k][1] - flux[i][j][k][1] ) + dx2 * tx1 * ( u[i-1][j][k][1] - 2.0 * u[i][j][k][1] + u[i+1][j][k][1] );
	  rsd[i][j][k][2] = rsd[i][j][k][2] + tx3 * C3 * C4 * ( flux[i+1][j][k][2] - flux[i][j][k][2] ) + dx3 * tx1 * ( u[i-1][j][k][2] - 2.0 * u[i][j][k][2] + u[i+1][j][k][2] );
	  rsd[i][j][k][3] = rsd[i][j][k][3] + tx3 * C3 * C4 * ( flux[i+1][j][k][3] - flux[i][j][k][3] ) + dx4 * tx1 * ( u[i-1][j][k][3] - 2.0 * u[i][j][k][3] + u[i+1][j][k][3] );
	  rsd[i][j][k][4] = rsd[i][j][k][4] + tx3 * C3 * C4 * ( flux[i+1][j][k][4] - flux[i][j][k][4] ) + dx5 * tx1 * ( u[i-1][j][k][4] - 2.0 * u[i][j][k][4] + u[i+1][j][k][4] );
	}
      }

      // X-direction diffusion boundary updates (i=1, 2)
      // These updates apply to fixed i indices for varying j, k, m.
      // Parallelizable across j and k. Updates for different i are independent.
      for (m = 0; m < 5; m++) {
	rsd[1][j][k][m] = rsd[1][j][k][m] - dssp * ( + 5.0 * u[1][j][k][m] - 4.0 * u[2][j][k][m] + u[3][j][k][m] );
	rsd[2][j][k][m] = rsd[2][j][k][m] - dssp * ( - 4.0 * u[1][j][k][m] + 6.0 * u[2][j][k][m] - 4.0 * u[3][j][k][m] + u[4][j][k][m] );
      }

      // X-direction diffusion inner updates (i=3 to nx-4)
      // Depends on u[i-2...i+2]. Updates rsd[i] for i=3 to nx-4. Parallelizable across i.
      ist1 = 3; iend1 = nx - 4;
      for (i = ist1; i <= iend1; i++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] = rsd[i][j][k][m] - dssp * ( u[i-2][j][k][m] - 4.0 * u[i-1][j][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i+1][j][k][m] + u[i+2][j][k][m] );
	}
      }

      // X-direction diffusion boundary updates (i=nx-3, nx-2)
      // These updates apply to fixed i indices for varying j, k, m.
      // Parallelizable across j and k. Updates for different i are independent.
      for (m = 0; m < 5; m++) {
	rsd[nx-3][j][k][m] = rsd[nx-3][j][k][m] - dssp * ( u[nx-5][j][k][m] - 4.0 * u[nx-4][j][k][m] + 6.0 * u[nx-3][j][k][m] - 4.0 * u[nx-2][j][k][m]  );
	rsd[nx-2][j][k][m] = rsd[nx-2][j][k][m] - dssp * ( u[nx-4][j][k][m] - 4.0 * u[nx-3][j][k][m] + 5.0 * u[nx-2][j][k][m] );
      }
    } // end for k
  } // end for j

  //---------------------------------------------------------------------
  // Phase 3: Y-direction calculations
  // Computes fluxes in the y-direction and updates rsd based on them.
  // Requires u. Overwrites flux. Updates rsd.
  // Implicit barrier required before this phase due to 'flux' reuse.
  // Inside this phase, many computations are parallelizable.
  //---------------------------------------------------------------------

  // Sub-phase 3.1: Calculate convective flux based on u (y-direction)
  // flux[i][j][k][m] depends only on u[i][j][k][m]. Overwrites flux from Phase 2.
  // Parallelizable across i, j, and k. Loop structure is i, j, k.
  // Data Dependencies: Reads u, Writes flux. No loop-carried dependencies.
  // Load Balance: Good if iterations are distributed evenly.
  // Synchronization: Not required within this sub-phase.
  // (L1 and L2 are used as in original code, effectively 0 to ny-1 for j)
  L1 = 0; L2 = ny-1;
  for (i = ist; i <= iend; i++) { // Note different i bounds compared to X-flux calc
    for (j = L1; j <= L2; j++) {
      for (k = 1; k <= nz - 2; k++) {
	flux[i][j][k][0] = u[i][j][k][2]; // Uses u from index 2
	u31 = u[i][j][k][2] / u[i][j][k][0]; // Uses u from index 2
	q = 0.50 * (  u[i][j][k][1] * u[i][j][k][1]
		      + u[i][j][k][2] * u[i][j][k][2]
		      + u[i][j][k][3] * u[i][j][k][3] )
	  / u[i][j][k][0];
	flux[i][j][k][1] = u[i][j][k][1] * u31;
	flux[i][j][k][2] = u[i][j][k][2] * u31 + C2 * (u[i][j][k][4]-q); // Uses u from index 2
	flux[i][j][k][3] = u[i][j][k][3] * u31;
	flux[i][j][k][4] = ( C1 * u[i][j][k][4] - C2 * q ) * u31;
      }
    }
  }

  // Sub-phase 3.2: Apply y-direction contributions to rsd
  // Updates rsd based on y-direction fluxes and diffusion.
  // Loops are structured with i, k outermost, then j, m inner.
  // Parallelizable across i and k for the outer loops, and across j for the inner loops.
  // Data Dependencies: Reads rsd, u, flux (calculated in Sub-phase 3.1). Writes rsd.
  // Updates to rsd[i][j][k][m] depend on flux[i][j-1/j+1][k][m] or flux[i][j/j+1][k][m]
  // and u[i][j-2..j+2][k][m]. No loop-carried dependencies on rsd or flux within the j loops.
  // Load Balance: Can be good depending on how loops (i, k, j) are parallelized.
  // Synchronization: Not required within this sub-phase if parallelized carefully. Implicit sync needed before Phase 4.
  for (i = ist; i <= iend; i++) {
    for (k = 1; k <= nz - 2; k++) { // Note original order k, j here vs j, k in X-dir rsd updates

      // rsd update based on convective flux differences in y-direction
      // Reads flux calculated in Sub-phase 3.1. Parallelizable across j.
      for (j = jst; j <= jend; j++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] =  rsd[i][j][k][m]
	    - ty2 * ( flux[i][j+1][k][m] - flux[i][j-1][k][m] );
	}
      }

      // Calculate viscous flux in y-direction (overwrites previous flux)
      // Depends on u[j] and u[j-1]. Parallelizable across j.
      // (L2 is used as in original code, effectively ny-1 for j)
      L2 = ny-1;
      for (j = jst; j <= L2; j++) {
	tmp = 1.0 / u[i][j][k][0];
	u21j = tmp * u[i][j][k][1]; u31j = tmp * u[i][j][k][2];
	u41j = tmp * u[i][j][k][3]; u51j = tmp * u[i][j][k][4];
	tmp = 1.0 / u[i][j-1][k][0];
	u21jm1 = tmp * u[i][j-1][k][1]; u31jm1 = tmp * u[i][j-1][k][2];
	u41jm1 = tmp * u[i][j-1][k][3]; u51jm1 = tmp * u[i][j-1][k][4];
	flux[i][j][k][1] = ty3 * ( u21j - u21jm1 );
	flux[i][j][k][2] = (4.0/3.0) * ty3 * (u31j-u31jm1); // Note: (4.0/3.0) on u31j here
	flux[i][j][k][3] = ty3 * ( u41j - u41jm1 );
	// Assuming pow2 is defined elsewhere (e.g., #define pow2(a) ((a)*(a)))
	flux[i][j][k][4] = 0.50 * ( 1.0 - C1*C5 ) * ty3 * ( ( pow2(u21j) + pow2(u31j) + pow2(u41j) ) - ( pow2(u21jm1) + pow2(u31jm1) + pow2(u41jm1) ) ) + (1.0/6.0) * ty3 * ( pow2(u31j) - pow2(u31jm1) ) + C1 * C5 * ty3 * ( u51j - u51jm1 ); // Note: pow2(u31j) here
      }

      // rsd update based on diffusion terms in y-direction and viscous flux
      // Uses flux calculated in the previous loop. Parallelizable across j.
      for (j = jst; j <= jend; j++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][0] = rsd[i][j][k][0] + dy1 * ty1 * ( u[i][j-1][k][0] - 2.0 * u[i][j][k][0] + u[i][j+1][k][0] );
	  rsd[i][j][k][1] = rsd[i][j][k][1] + ty3 * C3 * C4 * ( flux[i][j+1][k][1] - flux[i][j][k][1] ) + dy2 * ty1 * ( u[i][j-1][k][1] - 2.0 * u[i][j][k][1] + u[i][j+1][k][1] );
	  rsd[i][j][k][2] = rsd[i][j][k][2] + ty3 * C3 * C4 * ( flux[i][j+1][k][2] - flux[i][j][k][2] ) + dy3 * ty1 * ( u[i][j-1][k][2] - 2.0 * u[i][j][k][2] + u[i][j+1][k][2] );
	  rsd[i][j][k][3] = rsd[i][j][k][3] + ty3 * C3 * C4 * ( flux[i][j+1][k][3] - flux[i][j][k][3] ) + dy4 * ty1 * ( u[i][j-1][k][3] - 2.0 * u[i][j][k][3] + u[i][j+1][k][3] );
	  rsd[i][j][k][4] = rsd[i][j][k][4] + ty3 * C3 * C4 * ( flux[i][j+1][k][4] - flux[i][j][k][4] ) + dy5 * ty1 * ( u[i][j-1][k][4] - 2.0 * u[i][j][k][4] + u[i][j+1][k][4] );
	}
      }

      // Y-direction diffusion boundary updates (j=1, 2)
      // These updates apply to fixed j indices for varying i, k, m.
      // Parallelizable across i and k. Updates for different j are independent.
      for (m = 0; m < 5; m++) {
	rsd[i][1][k][m] = rsd[i][1][k][m] - dssp * ( + 5.0 * u[i][1][k][m] - 4.0 * u[i][2][k][m] + u[i][3][k][m] );
	rsd[i][2][k][m] = rsd[i][2][k][m] - dssp * ( - 4.0 * u[i][1][k][m] + 6.0 * u[i][2][k][m] - 4.0 * u[i][3][k][m] + u[i][4][k][m] );
      }

      // Y-direction diffusion inner updates (j=3 to ny-4)
      // Depends on u[j-2...j+2]. Updates rsd[j] for j=3 to ny-4. Parallelizable across j.
      jst1 = 3; jend1 = ny - 4;
      for (j = jst1; j <= jend1; j++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] = rsd[i][j][k][m] - dssp * ( u[i][j-2][k][m] - 4.0 * u[i][j-1][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j+1][k][m] + u[i][j+2][k][m] );
	}
      }

      // Y-direction diffusion boundary updates (j=ny-3, ny-2)
      // These updates apply to fixed j indices for varying i, k, m.
      // Parallelizable across i and k. Updates for different j are independent.
      for (m = 0; m < 5; m++) {
	rsd[i][ny-3][k][m] = rsd[i][ny-3][k][m] - dssp * ( u[i][ny-5][k][m] - 4.0 * u[i][ny-4][k][m] + 6.0 * u[i][ny-3][k][m] - 4.0 * u[i][ny-2][k][m]  );
	rsd[i][ny-2][k][m] = rsd[i][ny-2][k][m] - dssp * ( u[i][ny-4][k][m] - 4.0 * u[i][ny-3][k][m] + 5.0 * u[i][ny-2][k][m] );
      }
    } // end for k
  } // end for i

  //---------------------------------------------------------------------
  // Phase 4: Z-direction calculations
  // Computes fluxes in the z-direction and updates rsd based on them.
  // Requires u. Overwrites flux. Updates rsd.
  // Implicit barrier required before this phase due to 'flux' reuse.
  // Inside this phase, many computations are parallelizable.
  //---------------------------------------------------------------------

  // Sub-phase 4.1: Calculate convective flux based on u (z-direction)
  // flux[i][j][k][m] depends only on u[i][j][k][m]. Overwrites flux from Phase 3.
  // Parallelizable across i, j, and k. Loop structure is i, j, k.
  // Data Dependencies: Reads u, Writes flux. No loop-carried dependencies.
  // Load Balance: Good if iterations are distributed evenly.
  // Synchronization: Not required within this sub-phase.
  for (i = ist; i <= iend; i++) { // Note different i bounds compared to X-flux calc (0..nx-1)
    for (j = jst; j <= jend; j++) { // Note different j bounds compared to Y-flux calc (0..ny-1)
      for (k = 0; k <= nz-1; k++) { // Note k bounds are 0..nz-1
	flux[i][j][k][0] = u[i][j][k][3]; // Uses u from index 3
	u41 = u[i][j][k][3] / u[i][j][k][0]; // Uses u from index 3
	q = 0.50 * (  u[i][j][k][1] * u[i][j][k][1]
		      + u[i][j][k][2] * u[i][j][k][2]
		      + u[i][j][k][3] * u[i][j][k][3] )
	  / u[i][j][k][0];
	flux[i][j][k][1] = u[i][j][k][1] * u41;
	flux[i][j][k][2] = u[i][j][k][2] * u41;
	flux[i][j][k][3] = u[i][j][k][3] * u41 + C2 * (u[i][j][k][4]-q); // Uses u from index 3
	flux[i][j][k][4] = ( C1 * u[i][j][k][4] - C2 * q ) * u41;
      }

      // Sub-phase 4.2: Apply z-direction contributions to rsd
      // Updates rsd based on z-direction fluxes and diffusion.
      // Loops are structured with i, j outermost, then k, m inner.
      // Parallelizable across i and j for the outer loops, and across k for the inner loops.
      // Data Dependencies: Reads rsd, u, flux (calculated in Sub-phase 4.1). Writes rsd.
      // Updates to rsd[i][j][k][m] depend on flux[i][j][k-1/k+1][m] or flux[i][j][k/k+1][m]
      // and u[i][j][k-2..k+2][m]. No loop-carried dependencies on rsd or flux within the k loops.
      // Load Balance: Can be good depending on how loops (i, j, k) are parallelized.
      // Synchronization: Not required within this sub-phase if parallelized carefully. No implicit sync needed after this phase within rhs.

      // rsd update based on convective flux differences in z-direction
      // Reads flux calculated in Sub-phase 4.1. Parallelizable across k.
      // Original loop is k=1 to nz-2.
      for (k = 1; k <= nz - 2; k++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] =  rsd[i][j][k][m]
	    - tz2 * ( flux[i][j][k+1][m] - flux[i][j][k-1][m] );
	}
      }

      // Calculate viscous flux in z-direction (overwrites previous flux)
      // Depends on u[k] and u[k-1]. Parallelizable across k.
      // Original loop is k=1 to nz-1.
      for (k = 1; k <= nz-1; k++) {
	tmp = 1.0 / u[i][j][k][0];
	u21k = tmp * u[i][j][k][1]; u31k = tmp * u[i][j][k][2];
	u41k = tmp * u[i][j][k][3]; u51k = tmp * u[i][j][k][4];
	tmp = 1.0 / u[i][j][k-1][0];
	u21km1 = tmp * u[i][j][k-1][1]; u31km1 = tmp * u[i][j][k-1][2];
	u41km1 = tmp * u[i][j][k-1][3]; u51km1 = tmp * u[i][j][k-1][4];
	flux[i][j][k][1] = tz3 * ( u21k - u21km1 );
	flux[i][j][k][2] = tz3 * ( u31k - u31km1 );
	flux[i][j][k][3] = (4.0/3.0) * tz3 * (u41k-u41km1); // Note: (4.0/3.0) on u41k here
	// Assuming pow2 is defined elsewhere (e.g., #define pow2(a) ((a)*(a)))
	flux[i][j][k][4] = 0.50 * ( 1.0 - C1*C5 ) * tz3 * ( ( pow2(u21k) + pow2(u31k) + pow2(u41k) ) - ( pow2(u21km1) + pow2(u31km1) + pow2(u41km1) ) ) + (1.0/6.0) * tz3 * ( pow2(u41k) - pow2(u41km1) ) + C1 * C5 * tz3 * ( u51k - u51km1 ); // Note: pow2(u41k) here
      }

      // rsd update based on diffusion terms in z-direction and viscous flux
      // Uses flux calculated in the previous loop. Parallelizable across k.
      // Original loop is k=1 to nz-2.
      for (k = 1; k <= nz - 2; k++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][0] = rsd[i][j][k][0] + dz1 * tz1 * ( u[i][j][k-1][0] - 2.0 * u[i][j][k][0] + u[i][j][k+1][0] );
	  rsd[i][j][k][1] = rsd[i][j][k][1] + tz3 * C3 * C4 * ( flux[i][j][k+1][1] - flux[i][j][k][1] ) + dz2 * tz1 * ( u[i][j][k-1][1] - 2.0 * u[i][j][k][1] + u[i][j][k+1][1] );
	  rsd[i][j][k][2] = rsd[i][j][k][2] + tz3 * C3 * C4 * ( flux[i][j][k+1][2] - flux[i][j][k][2] ) + dz3 * tz1 * ( u[i][j][k-1][2] - 2.0 * u[i][j][k][2] + u[i][j][k+1][2] );
	  rsd[i][j][k][3] = rsd[i][j][k][3] + tz3 * C3 * C4 * ( flux[i][j][k+1][3] - flux[i][j][k][3] ) + dz4 * tz1 * ( u[i][j][k-1][3] - 2.0 * u[i][j][k][3] + u[i][j][k+1][3] );
	  rsd[i][j][k][4] = rsd[i][j][k][4] + tz3 * C3 * C4 * ( flux[i][j][k+1][4] - flux[i][j][k][4] ) + dz5 * tz1 * ( u[i][j][k-1][4] - 2.0 * u[i][j][k][4] + u[i][j][k+1][4] );
	}
      }

      // Z-direction diffusion boundary updates (k=1, 2)
      // These updates apply to fixed k indices for varying i, j, m.
      // Parallelizable across i and j. Updates for different k are independent.
      for (m = 0; m < 5; m++) {
	rsd[i][j][1][m] = rsd[i][j][1][m] - dssp * ( + 5.0 * u[i][j][1][m] - 4.0 * u[i][j][2][m] + u[i][j][3][m] );
	rsd[i][j][2][m] = rsd[i][j][2][m] - dssp * ( - 4.0 * u[i][j][1][m] + 6.0 * u[i][j][2][m] - 4.0 * u[i][j][3][m] + u[i][j][4][m] );
      }

      // Z-direction diffusion inner updates (k=3 to nz-4)
      // Depends on u[k-2...k+2]. Updates rsd[k] for k=3 to nz-4. Parallelizable across k.
      // Original loop is k=3 to nz-4.
      for (k = 3; k <= nz - 4; k++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] = rsd[i][j][k][m] - dssp * ( u[i][j][k-2][m] - 4.0 * u[i][j][k-1][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j][k+1][m] + u[i][j][k+2][m] );
	}
      }

      // Z-direction diffusion boundary updates (k=nz-3, nz-2)
      // These updates apply to fixed k indices for varying i, j, m.
      // Parallelizable across i and j. Updates for different k are independent.
      for (m = 0; m < 5; m++) {
	rsd[i][j][nz-3][m] = rsd[i][j][nz-3][m] - dssp * ( u[i][j][nz-5][m] - 4.0 * u[i][j][nz-4][m] + 6.0 * u[i][j][nz-3][m] - 4.0 * u[i][j][nz-2][m]  );
	rsd[i][j][nz-2][m] = rsd[i][j][nz-2][m] - dssp * ( u[i][j][nz-4][m] - 4.0 * u[i][j][nz-3][m] + 5.0 * u[i][j][nz-2][m] );
      }
    } // end for j
  } // end for i
}