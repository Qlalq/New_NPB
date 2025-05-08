static void verify(int no_time_steps, char *class, boolean *verified) {
  double xcrref[5],xceref[5],xcrdif[5],xcedif[5],
    epsilon, xce[5], xcr[5], dtref; // Declare variables at start

  // --- Section 1: Initial Computations and Setup ---
  // This section performs initial calculations and sets up scalar values.
  // It includes calls to external functions, which might contain internal parallelism,
  // but from the perspective of 'verify', these are sequential dependencies.
  epsilon = 1.0e-08;

  error_norm(xce); // Assume populates xce
  compute_rhs();   // Assume prepares data for rhs_norm
  rhs_norm(xcr);   // Assume populates xcr

  // --- Section 2: Parallelizable Array Computations (Fixed size 5 loops) ---
  // These loops iterate over a small, fixed-size array (size 5).
  // Each iteration is independent, making them candidates for OpenMP parallel for
  // if the array size was significantly larger.

  // Loop 1: Scale xcr array elements. Iterations are independent.
  int m; // Declare loop variable m here
  for (m = 0; m < 5; m++) {
    xcr[m] = xcr[m] / dt; // dt is a scalar or constant here
  }

  // Initialize verification status. These are scalar writes.
  *class = 'U';
  *verified = TRUE;

  // Loop 2: Initialize reference arrays. Iterations are independent.
  // This initialization is overwritten later if a specific class is matched.
  for (m = 0; m < 5; m++) {
    xcrref[m] = 1.0;
    xceref[m] = 1.0;
  }

  // --- Section 3: Serial Configuration based on Input Parameters ---
  // This section determines the specific reference values based on problem size
  // and time steps. It's a control flow structure and is inherently serial.
  if ( grid_points[0] == 12 &&
       grid_points[1] == 12 &&
       grid_points[2] == 12 &&
       no_time_steps == 100) {
    *class = 'S';
    dtref = 1.5e-2;
    xcrref[0] = 2.7470315451339479e-02;
    xcrref[1] = 1.0360746705285417e-02;
    xcrref[2] = 1.6235745065095532e-02;
    xcrref[3] = 1.5840557224455615e-02;
    xcrref[4] = 3.4849040609362460e-02;
    xceref[0] = 2.7289258557377227e-05;
    xceref[1] = 1.0364446640837285e-05;
    xceref[2] = 1.6154798287166471e-05;
    xceref[3] = 1.5750704994480102e-05;
    xceref[4] = 3.4177666183390531e-05;
  } else if (grid_points[0] == 36 &&
	     grid_points[1] == 36 &&
	     grid_points[2] == 36 &&
	     no_time_steps == 400) {
    *class = 'W';
    dtref = 1.5e-3;
    xcrref[0] = 0.1893253733584e-02;
    xcrref[1] = 0.1717075447775e-03;
    xcrref[2] = 0.2778153350936e-03;
    xcrref[3] = 0.2887475409984e-03;
    xcrref[4] = 0.3143611161242e-02;
    xceref[0] = 0.7542088599534e-04;
    xceref[1] = 0.6512852253086e-05;
    xceref[2] = 0.1049092285688e-04;
    xceref[3] = 0.1128838671535e-04;
    xceref[4] = 0.1212845639773e-03;
  } else if (grid_points[0] == 64 &&
	     grid_points[1] == 64 &&
	     grid_points[2] == 64 &&
	     no_time_steps == 400 ) {
    *class = 'A';
    dtref = 1.5e-3;
    xcrref[0] = 2.4799822399300195;
    xcrref[1] = 1.1276337964368832;
    xcrref[2] = 1.5028977888770491;
    xcrref[3] = 1.4217816211695179;
    xcrref[4] = 2.1292113035138280;
    xceref[0] = 1.0900140297820550e-04;
    xceref[1] = 3.7343951769282091e-05;
    xceref[2] = 5.0092785406541633e-05;
    xceref[3] = 4.7671093939528255e-05;
    xceref[4] = 1.3621613399213001e-04;
  } else if (grid_points[0] == 102 &&
	     grid_points[1] == 102 &&
	     grid_points[2] == 102 &&
	     no_time_steps == 400) {
    *class = 'B';
    dtref = 1.0e-3;
    xcrref[0] = 0.6903293579998e+02;
    xcrref[1] = 0.3095134488084e+02;
    xcrref[2] = 0.4103336647017e+02;
    xcrref[3] = 0.3864769009604e+02;
    xcrref[4] = 0.5643482272596e+02;
    xceref[0] = 0.9810006190188e-02;
    xceref[1] = 0.1022827905670e-02;
    xceref[2] = 0.1720597911692e-02;
    xceref[3] = 0.1694479428231e-02;
    xceref[4] = 0.1847456263981e-01;
  } else if (grid_points[0] == 162 &&
	     grid_points[1] == 162 &&
	     grid_points[2] == 162 &&
	     no_time_steps == 400) {
    *class = 'C';
    dtref = 0.67e-3;
    xcrref[0] = 0.5881691581829e+03;
    xcrref[1] = 0.2454417603569e+03;
    xcrref[2] = 0.3293829191851e+03;
    xcrref[3] = 0.3081924971891e+03;
    xcrref[4] = 0.4597223799176e+03;
    xceref[0] = 0.2598120500183e+00;
    xceref[1] = 0.2590888922315e-01;
    xceref[2] = 0.5132886416320e-01;
    xceref[3] = 0.4806073419454e-01;
    xceref[4] = 0.5483377491301e+00;
  } else {
    // If class remains 'U' because no size matched, verification is not possible
    // against reference values. Set verified to FALSE explicitly here.
    *verified = FALSE;
  }

  // --- Section 4: Parallelizable Array Computation (Fixed size 5 loop) ---
  // This loop computes difference ratios for verification.
  // Each iteration is independent, making it a candidate for OpenMP parallel for
  // if the array size was significantly larger.

  // Loop 3: Compute difference ratios. Iterations are independent.
  for (m = 0; m < 5; m++) {
    // These computations depend on xcr, xce (from Section 1 & 2)
    // and xcrref, xceref (from Section 2 & 3).
    xcrdif[m] = fabs((xcr[m]-xcrref[m])/xcrref[m]) ;
    xcedif[m] = fabs((xce[m]-xceref[m])/xceref[m]);
  }

  // --- Section 5: Serial Verification and Output ---
  // This section performs the actual verification checks and prints results.
  // It involves conditional logic, I/O operations (printf), and updating a
  // shared variable (*verified). I/O is a synchronization point, and updates
  // to *verified would require careful handling (e.g., critical section or atomic)
  // if this section were parallelized. It's generally simpler and safer to keep
  // this part serial after all computations are complete.

  if (*class != 'U') {
    printf(" Verification being performed for class %1c\n", *class);
    printf(" accuracy setting for epsilon = %20.13e\n", epsilon);
    // Check dt against reference if class is known.
    if (fabs(dt-dtref) > epsilon) {
      *verified = FALSE; // Update shared flag on failure
      *class = 'U';      // Indicate verification failed by setting class to 'U'
      printf(" DT does not match the reference value of %15.8e\n", dtref);
    }
  } else {
    printf(" Unknown class\n");
  }

  // Print header for xcr comparison results.
  if (*class != 'U') {
    printf(" Comparison of RMS-norms of residual\n");
  } else {
    printf(" RMS-norms of residual\n");
  }

  // Loop 4: Print xcr comparison results and update verification flag.
  // Contains I/O and potential shared write. Best kept serial.
  for (m = 0; m < 5; m++) {
    if (*class == 'U') {
      printf("          %2d%20.13e\n", m, xcr[m]);
    } else if (xcrdif[m] > epsilon) {
      *verified = FALSE; // Update shared flag on failure
      printf(" FAILURE: %2d%20.13e%20.13e%20.13e\n",
	     m,xcr[m],xcrref[m],xcrdif[m]);
    } else {
      printf("          %2d%20.13e%20.13e%20.13e\n",
	     m,xcr[m],xcrref[m],xcrdif[m]);
    }
  }

  // Print header for xce comparison results.
  if (*class != 'U') {
    printf(" Comparison of RMS-norms of solution error\n");
  } else {
    printf(" RMS-norms of solution error\n");
  }

  // Loop 5: Print xce comparison results and update verification flag.
  // Contains I/O and potential shared write. Best kept serial.
  for (m = 0; m < 5; m++) {
    if (*class == 'U') {
      printf("          %2d%20.13e\n", m, xce[m]);
    } else if (xcedif[m] > epsilon) {
      *verified = FALSE; // Update shared flag on failure
      printf(" FAILURE: %2d%20.13e%20.13e%20.13e\n",
	     m,xce[m],xceref[m],xcedif[m]);
    } else {
      printf("          %2d%20.13e%20.13e%20.13e\n",
	     m,xce[m],xceref[m],xcedif[m]);
    }
  }

  // Print final verification summary.
  if (*class == 'U') {
    printf(" No reference values provided\n");
    printf(" No verification performed\n");
  } else if (*verified) {
    printf(" Verification Successful\n");
  } else {
    printf(" Verification failed\n");
  }
}