// Note: This code structure is suitable for OpenMP parallelization.
// The following comments highlight parallelizable loops and data considerations
// based on the analysis of data dependencies, load balance, and synchronization points.

// Assume ue, buf, cuf, q are declared in the calling scope with sufficient size
// to hold data for any dimension (grid_points[0], grid_points[1], or grid_points[2]).
// e.g., double ue[5][MAX_GRID_DIMENSION], buf[5][MAX_GRID_DIMENSION], cuf[MAX_GRID_DIMENSION], q[MAX_GRID_DIMENSION];
// These temporary arrays need to be thread-private when parallelizing the loops below.

// Assume grid_points and constants like dnzm1, dnym1, dnxm1, tx2, etc.,
// dssp, c1, c2 are accessible (globals or passed in).
// Assume exact_solution function is defined and accessible.

static void exact_rhs(void) {
  // Variables used within loops. These are typically thread-private
  // when the loops they are in are parallelized.
  double dtemp[5], xi, eta, zeta, dtpp;
  int m, i, j, k, ip1, im1, jp1, jm1, km1, kp1;

  // Section 1: Initialization of forcing array
  // This section is highly parallelizable. Each iteration writes to a unique element
  // of the forcing array. No dependencies between iterations.
  // Potential OpenMP parallelization: distribute iterations across threads.
  // Example directive: #pragma omp parallel for collapse(4) private(m,i,j,k)
  for (m = 0; m < 5; m++) {
    for (i = 0; i <= grid_points[0]-1; i++) {
      for (j = 0; j <= grid_points[1]-1; j++) {
	for (k= 0; k <= grid_points[2]-1; k++) {
	  forcing[m][i][j][k] = 0.0;
	}
      }
    }
  }

  // Section 2: X-direction computations (updates based on i-stencil)
  // Parallelization strategy: Parallelize the outer loops over k and j.
  // For each (j, k) pair, the inner loops iterate over 'i'.
  // The computations within the inner 'i' loops for a given (j, k) use temporary
  // 1D arrays (ue, buf, cuf, q) indexed by 'i'.
  // These temporary arrays (sized grid_points[0]) must be private to each thread
  // executing a (j, k) iteration block to avoid data races.
  // Updates to forcing[m][i][j][k] are independent for different (i, j, k) tuples,
  // so no synchronization is needed for writing to 'forcing' within this section
  // if parallelized over j and/or k.
  // Potential OpenMP parallelization: distribute (j,k) iterations across threads.
  // Example directive: #pragma omp parallel for collapse(2) private(j,k, zeta, eta, xi, dtemp, dtpp, im1, ip1, m, ue, buf, cuf, q)
  for (k = 1; k <= grid_points[2]-2; k++) {
    zeta = (double)k * dnzm1; // Private variable
    for (j = 1; j <= grid_points[1]-2; j++) {
      eta = (double)j * dnym1; // Private variable

      // These loops process the current (j, k) plane slice of the grid.
      // Temporary arrays (ue, buf, cuf, q indexed by i) are specific to this (j, k).
      // They must be private to the thread executing this block.

      // Loop 1 (i): Compute temporary values (ue, buf, cuf, q) indexed by i.
      for (i = 0; i <= grid_points[0]-1; i++) {
	xi = (double)i * dnxm1; // Private variable
	exact_solution(xi, eta, zeta, dtemp); // dtemp is private/local
	for (m = 0; m < 5; m++) { // m is private
	  ue[m][i] = dtemp[m]; // Writes to private ue[m][i]
	}
	dtpp = 1.0 / dtemp[0]; // Private variable
	for (m = 1; m < 5; m++) { // m is private
	  buf[m][i] = dtpp * dtemp[m]; // Writes to private buf[m][i]
	}
	cuf[i] = buf[1][i] * buf[1][i]; // Writes to private cuf[i]
	buf[0][i] = cuf[i] + buf[2][i] * buf[2][i] + buf[3][i] * buf[3][i]; // Writes to private buf[0][i]
	q[i] = 0.5 * (buf[1][i]*ue[1][i] + buf[2][i]*ue[2][i]
		      + buf[3][i]*ue[3][i]); // Writes to private q[i]
      }

      // Loop 2 (i): Update forcing based on temporary values using an i-stencil.
      // Reads private ue/buf/cuf/q[i-1..i+1], writes to forcing[m][i][j][k].
      // The write to forcing[m][i][j][k] is unique for this (i,j,k).
      for (i = 1; i <= grid_points[0]-2; i++) {
	im1 = i-1; // Private variable
	ip1 = i+1; // Private variable
	forcing[0][i][j][k] = forcing[0][i][j][k] -
	  tx2*( ue[1][ip1]-ue[1][im1] )+
	  dx1tx1*(ue[0][ip1]-2.0*ue[0][i]+ue[0][im1]);
	forcing[1][i][j][k] = forcing[1][i][j][k]
	  - tx2 * ((ue[1][ip1]*buf[1][ip1]+c2*(ue[4][ip1]-q[ip1]))-
                   (ue[1][im1]*buf[1][im1]+c2*(ue[4][im1]-q[im1])))+
	  xxcon1*(buf[1][ip1]-2.0*buf[1][i]+buf[1][im1])+
	  dx2tx1*( ue[1][ip1]-2.0* ue[1][i]+ue[1][im1]);
	forcing[2][i][j][k] = forcing[2][i][j][k]
	  - tx2 * (ue[2][ip1]*buf[1][ip1]-ue[2][im1]*buf[1][im1])+
	  xxcon2*(buf[2][ip1]-2.0*buf[2][i]+buf[2][im1])+
	  dx3tx1*( ue[2][ip1]-2.0*ue[2][i] +ue[2][im1]);
	forcing[3][i][j][k] = forcing[3][i][j][k]
	  - tx2*(ue[3][ip1]*buf[1][ip1]-ue[3][im1]*buf[1][im1])+
	  xxcon2*(buf[3][ip1]-2.0*buf[3][i]+buf[3][im1])+
	  dx4tx1*( ue[3][ip1]-2.0* ue[3][i]+ ue[3][im1]);
	forcing[4][i][j][k] = forcing[4][i][j][k]
	  - tx2*(buf[1][ip1]*(c1*ue[4][ip1]-c2*q[ip1])-
		 buf[1][im1]*(c1*ue[4][im1]-c2*q[im1]))+
	  0.5*xxcon3*(buf[0][ip1]-2.0*buf[0][i]+
		      buf[0][im1])+
	  xxcon4*(cuf[ip1]-2.0*cuf[i]+cuf[im1])+
	  xxcon5*(buf[4][ip1]-2.0*buf[4][i]+buf[4][im1])+
	  dx5tx1*( ue[4][ip1]-2.0* ue[4][i]+ ue[4][im1]);
      }
      // Boundary updates (i=1, 2, gp[0]-3, gp[0]-2 and range 3 to gp[0]-4)
      // These also read private ue/buf/cuf/q and write to unique forcing elements.
      for (m = 0; m < 5; m++) { // m is private
	i = 1; // i is private here
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (5.0*ue[m][i] - 4.0*ue[m][i+1] +ue[m][i+2]);
	i = 2; // i is private here
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (-4.0*ue[m][i-1] + 6.0*ue[m][i] -
 	    4.0*ue[m][i+1] +     ue[m][i+2]);
      }
      for (m = 0; m < 5; m++) { // m is private
	for (i = 3; i <= grid_points[0]-4; i++) { // i is private
	  forcing[m][i][j][k] = forcing[m][i][j][k] - dssp*
	    (ue[m][i-2] - 4.0*ue[m][i-1] +
	     6.0*ue[m][i] - 4.0*ue[m][i+1] + ue[m][i+2]);
	}
      }
      for (m = 0; m < 5; m++) { // m is private
	i = grid_points[0]-3; // i is private here
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (ue[m][i-2] - 4.0*ue[m][i-1] +
	   6.0*ue[m][i] - 4.0*ue[m][i+1]);
	i = grid_points[0]-2; // i is private here
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (ue[m][i-2] - 4.0*ue[m][i-1] + 5.0*ue[m][i]);
      }
    } // End of loop j
  } // End of loop k

  // Section 3: Y-direction computations (updates based on j-stencil)
  // Parallelization strategy: Parallelize the outer loops over k and i.
  // For each (i, k) pair, the inner loops iterate over 'j'.
  // Temporary arrays ue, buf, cuf, q (indexed by j) must be private to each thread.
  // Variables zeta, xi, eta, dtemp, dtpp, jm1, jp1, m are temporary/loop private.
  // Updates to forcing[m][i][j][k] are independent for different (i, j, k) tuples,
  // so no synchronization is needed for writing to 'forcing' within this section
  // if parallelized over i and/or k.
  // Potential OpenMP parallelization: distribute (i,k) iterations across threads.
  // Example directive: #pragma omp parallel for collapse(2) private(i,k, zeta, xi, eta, dtemp, dtpp, jm1, jp1, m, ue, buf, cuf, q)
  for (k = 1; k <= grid_points[2]-2; k++) {
    zeta = (double)k * dnzm1; // Private variable
    for (i = 1; i <= grid_points[0]-2; i++) {
      xi = (double)i * dnxm1; // Private variable

      // These loops process the current (i, k) plane slice of the grid.
      // Temporary arrays (ue, buf, cuf, q indexed by j) are specific to this (i, k).
      // They must be private to the thread executing this block.

      // Loop 1 (j): Compute temporary values (ue, buf, cuf, q) indexed by j.
      for (j = 0; j <= grid_points[1]-1; j++) {
	eta = (double)j * dnym1; // Private variable
	exact_solution(xi, eta, zeta, dtemp); // dtemp is private/local
	for (m = 0; m < 5; m++) { // m is private
	  ue[m][j] = dtemp[m]; // Writes to private ue[m][j]
	}
	dtpp = 1.0/dtemp[0]; // Private variable
	for (m = 1; m < 5; m++) { // m is private
	  buf[m][j] = dtpp * dtemp[m]; // Writes to private buf[m][j]
	}
	cuf[j]   = buf[2][j] * buf[2][j]; // Writes to private cuf[j]
	buf[0][j] = cuf[j] + buf[1][j] * buf[1][j] +
	  buf[3][j] * buf[3][j]; // Writes to private buf[0][j]
	q[j] = 0.5*(buf[1][j]*ue[1][j] + buf[2][j]*ue[2][j] +
		    buf[3][j]*ue[3][j]); // Writes to private q[j]
      }

      // Loop 2 (j): Update forcing based on temporary values using a j-stencil.
      // Reads private ue/buf/cuf/q[j-1..j+1], writes to forcing[m][i][j][k].
      // The write to forcing[m][i][j][k] is unique for this (i,j,k).
      for (j = 1; j <= grid_points[1]-2; j++) {
	jm1 = j-1; // Private variable
	jp1 = j+1; // Private variable
	forcing[0][i][j][k] = forcing[0][i][j][k] -
	  ty2*( ue[2][jp1]-ue[2][jm1] )+
	  dy1ty1*(ue[0][jp1]-2.0*ue[0][j]+ue[0][jm1]);
	forcing[1][i][j][k] = forcing[1][i][j][k]
	  - ty2*(ue[1][jp1]*buf[2][jp1]-ue[1][jm1]*buf[2][jm1])+
	  yycon2*(buf[1][jp1]-2.0*buf[1][j]+buf[1][jm1])+
	  dy2ty1*( ue[1][jp1]-2.0* ue[1][j]+ ue[1][jm1]);
	forcing[2][i][j][k] = forcing[2][i][j][k]
	  - ty2*((ue[2][jp1]*buf[2][jp1]+c2*(ue[4][jp1]-q[jp1]))-
		 (ue[2][jm1]*buf[2][jm1]+c2*(ue[4][jm1]-q[jm1])))+
	  yycon1*(buf[2][jp1]-2.0*buf[2][j]+buf[2][jm1])+
	  dy3ty1*( ue[2][jp1]-2.0*ue[2][j] +ue[2][jm1]);
	forcing[3][i][j][k] = forcing[3][i][j][k]
	  - ty2*(ue[3][jp1]*buf[2][jp1]-ue[3][jm1]*buf[2][jm1])+
	  yycon2*(buf[3][jp1]-2.0*buf[3][j]+buf[3][jm1])+
	  dy4ty1*( ue[3][jp1]-2.0*ue[3][j]+ ue[3][jm1]);
	forcing[4][i][j][k] = forcing[4][i][j][k]
	  - ty2*(buf[2][jp1]*(c1*ue[4][jp1]-c2*q[jp1])-
		 buf[2][jm1]*(c1*ue[4][jm1]-c2*q[jm1]))+
	  0.5*yycon3*(buf[0][jp1]-2.0*buf[0][j]+
		      buf[0][jm1])+
	  yycon4*(cuf[jp1]-2.0*cuf[j]+cuf[jm1])+
	  yycon5*(buf[4][jp1]-2.0*buf[4][j]+buf[4][jm1])+
	  dy5ty1*(ue[4][jp1]-2.0*ue[4][j]+ue[4][jm1]);
      }
      // Boundary updates (j=1, 2, gp[1]-3, gp[1]-2 and range 3 to gp[1]-4)
      // These also read private ue/buf/cuf/q and write to unique forcing elements.
      for (m = 0; m < 5; m++) { // m is private
	j = 1; // j is private here
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (5.0*ue[m][j] - 4.0*ue[m][j+1] +ue[m][j+2]);
	j = 2; // j is private here
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (-4.0*ue[m][j-1] + 6.0*ue[m][j] -
	   4.0*ue[m][j+1] +       ue[m][j+2]);
      }
      for (m = 0; m < 5; m++) { // m is private
	for (j = 3; j <= grid_points[1]-4; j++) { // j is private
	  forcing[m][i][j][k] = forcing[m][i][j][k] - dssp*
	    (ue[m][j-2] - 4.0*ue[m][j-1] +
	     6.0*ue[m][j] - 4.0*ue[m][j+1] + ue[m][j+2]);
	}
      }
      for (m = 0; m < 5; m++) { // m is private
	j = grid_points[1]-3; // j is private here
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (ue[m][j-2] - 4.0*ue[m][j-1] +
	   6.0*ue[m][j] - 4.0*ue[m][j+1]);
	j = grid_points[1]-2; // j is private here
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (ue[m][j-2] - 4.0*ue[m][j-1] + 5.0*ue[m][j]);
      }
    } // End of loop i
  } // End of loop k

  // Section 4: Z-direction computations (updates based on k-stencil)
  // Parallelization strategy: Parallelize the outer loops over j and i.
  // For each (i, j) pair, the inner loops iterate over 'k'.
  // Temporary arrays ue, buf, cuf, q (indexed by k) must be private to each thread.
  // Variables eta, xi, zeta, dtemp, dtpp, km1, kp1, m are temporary/loop private.
  // Updates to forcing[m][i][j][k] are independent for different (i, j, k) tuples,
  // so no synchronization is needed for writing to 'forcing' within this section
  // if parallelized over i and/or j.
  // Potential OpenMP parallelization: distribute (i,j) iterations across threads.
  // Example directive: #pragma omp parallel for collapse(2) private(i,j, eta, xi, zeta, dtemp, dtpp, km1, kp1, m, ue, buf, cuf, q)
  for (j = 1; j <= grid_points[1]-2; j++) {
    eta = (double)j * dnym1; // Private variable
    for (i = 1; i <= grid_points[0]-2; i++) {
      xi = (double)i * dnxm1; // Private variable

      // These loops process the current (i, j) plane slice of the grid.
      // Temporary arrays (ue, buf, cuf, q indexed by k) are specific to this (i, j).
      // They must be private to the thread executing this block.

      // Loop 1 (k): Compute temporary values (ue, buf, cuf, q) indexed by k.
      for (k = 0; k <= grid_points[2]-1; k++) {
	zeta = (double)k * dnzm1; // Private variable
	exact_solution(xi, eta, zeta, dtemp); // dtemp is private/local
	for (m = 0; m < 5; m++) { // m is private
	  ue[m][k] = dtemp[m]; // Writes to private ue[m][k]
	}
	dtpp = 1.0/dtemp[0]; // Private variable
	for (m = 1; m < 5; m++) { // m is private
	  buf[m][k] = dtpp * dtemp[m]; // Writes to private buf[m][k]
	}
	cuf[k] = buf[3][k] * buf[3][k]; // Writes to private cuf[k]
	buf[0][k] = cuf[k] + buf[1][k] * buf[1][k] +
	  buf[2][k] * buf[2][k]; // Writes to private buf[0][k]
	q[k] = 0.5*(buf[1][k]*ue[1][k] + buf[2][k]*ue[2][k] +
		    buf[3][k]*ue[3][k]); // Writes to private q[k]
      }

      // Loop 2 (k): Update forcing based on temporary values using a k-stencil.
      // Reads private ue/buf/cuf/q[k-1..k+1], writes to forcing[m][i][j][k].
      // The write to forcing[m][i][j][k] is unique for this (i,j,k).
      for (k = 1; k <= grid_points[2]-2; k++) {
	km1 = k-1; // Private variable
	kp1 = k+1; // Private variable
	forcing[0][i][j][k] = forcing[0][i][j][k] -
	  tz2*( ue[3][kp1]-ue[3][km1] )+
	  dz1tz1*(ue[0][kp1]-2.0*ue[0][k]+ue[0][km1]);
	forcing[1][i][j][k] = forcing[1][i][j][k]
	  - tz2 * (ue[1][kp1]*buf[3][kp1]-ue[1][km1]*buf[3][km1])+
	  zzcon2*(buf[1][kp1]-2.0*buf[1][k]+buf[1][km1])+
	  dz2tz1*( ue[1][kp1]-2.0* ue[1][k]+ ue[1][km1]);
	forcing[2][i][j][k] = forcing[2][i][j][k]
	  - tz2 * (ue[2][kp1]*buf[3][kp1]-ue[2][km1]*buf[3][km1])+
	  zzcon2*(buf[2][kp1]-2.0*buf[2][k]+buf[2][km1])+
	  dz3tz1*(ue[2][kp1]-2.0*ue[2][k]+ue[2][km1]);
	forcing[3][i][j][k] = forcing[3][i][j][k]
	  - tz2 * ((ue[3][kp1]*buf[3][kp1]+c2*(ue[4][kp1]-q[kp1]))-
		   (ue[3][km1]*buf[3][km1]+c2*(ue[4][km1]-q[km1])))+
	  zzcon1*(buf[3][kp1]-2.0*buf[3][k]+buf[3][km1])+
	  dz4tz1*( ue[3][kp1]-2.0*ue[3][k] +ue[3][km1]);
	forcing[4][i][j][k] = forcing[4][i][j][k]
	  - tz2 * (buf[3][kp1]*(c1*ue[4][kp1]-c2*q[kp1])-
		   buf[3][km1]*(c1*ue[4][km1]-c2*q[km1]))+
	  0.5*zzcon3*(buf[0][kp1]-2.0*buf[0][k]
		      +buf[0][km1])+
	  zzcon4*(cuf[kp1]-2.0*cuf[k]+cuf[km1])+
	  zzcon5*(buf[4][kp1]-2.0*buf[4][k]+buf[4][km1])+
	  dz5tz1*( ue[4][kp1]-2.0* ue[4][k]+ ue[4][km1]);
      }
      // Boundary updates (k=1, 2, gp[2]-3, gp[2]-2 and range 3 to gp[2]-4)
      // These also read private ue/buf/cuf/q and write to unique forcing elements.
      for (m = 0; m < 5; m++) { // m is private
	k = 1; // k is private here
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (5.0*ue[m][k] - 4.0*ue[m][k+1] +ue[m][k+2]);
	k = 2; // k is private here
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (-4.0*ue[m][k-1] + 6.0*ue[m][k] -
	   4.0*ue[m][k+1] +       ue[m][k+2]);
      }
      for (m = 0; m < 5; m++) { // m is private
	for (k = 3; k <= grid_points[2]-4; k++) { // k is private
	  forcing[m][i][j][k] = forcing[m][i][j][k] - dssp*
	    (ue[m][k-2] - 4.0*ue[m][k-1] +
	     6.0*ue[m][k] - 4.0*ue[m][k+1] + ue[m][k+2]);
	}
      }
      for (m = 0; m < 5; m++) { // m is private
	k = grid_points[2]-3; // k is private here
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (ue[m][k-2] - 4.0*ue[m][k-1] +
	   6.0*ue[m][k] - 4.0*ue[m][k+1]);
	k = grid_points[2]-2; // k is private here
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (ue[m][k-2] - 4.0*ue[m][k-1] + 5.0*ue[m][k]);
      }
    } // End of loop i
  } // End of loop j

  // Section 5: Final sign change
  // This section is highly parallelizable. Each iteration writes to a unique element
  // of the forcing array. No dependencies between iterations.
  // Potential OpenMP parallelization: distribute iterations across threads.
  // Example directive: #pragma omp parallel for collapse(4) private(m,i,j,k)
  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = 1; k <= grid_points[2]-2; k++) {
	  forcing[m][i][j][k] = -1.0 * forcing[m][i][j][k];
	}
      }
    }
  }
}