#include <stdio.h>
#include <math.h> // Required for fabs

// Assume these are defined elsewhere (e.g., global variables or extern)
extern double dt;
extern int grid_points[3];
extern int TRUE; // Assume boolean TRUE is defined (e.g., #define TRUE 1)
extern int FALSE; // Assume boolean FALSE is defined (e.g., #define FALSE 0)
typedef int boolean; // Assume boolean type is defined

// Assume these functions are defined elsewhere
extern void error_norm(double *xce);
extern void compute_rhs();
extern void rhs_norm(double *xcr);

static void verify(int no_time_steps, char *class, boolean *verified) {
  double xcrref[5], xceref[5], xcrdif[5], xcedif[5],
         epsilon, xce[5], xcr[5], dtref;
  int m;

  // Phase 1: Initial setup, external calls, xcr scaling (Sequential)
  // This phase depends on external state/functions (error_norm, compute_rhs, rhs_norm, dt).
  // The loop scaling xcr is element-wise independent but small.
  epsilon = 1.0e-08;
  error_norm(xce); // External call, modifies xce
  compute_rhs();   // External call, modifies external state
  rhs_norm(xcr);     // External call, modifies xcr

  // Scale xcr array - Loop 1
  // Iterations are independent, operates on local xcr and global dt.
  for (m = 0; m < 5; m++) {
    xcr[m] = xcr[m] / dt;
  }

  // Phase 2: Set references and class (Sequential)
  // This phase reads global grid_points and no_time_steps, and writes
  // to pointer arguments (*class, *verified) and local reference arrays.
  // It determines the expected values for verification.
  *class = 'U'; // Initial assumption: Unknown class
  *verified = TRUE; // Initial assumption: Verified unless proven otherwise

  // Initialize reference arrays - Loop 2
  // Iterations are independent, operates on local xcrref, xceref.
  for (m = 0; m < 5; m++) {
    xcrref[m] = 1.0;
    xceref[m] = 1.0;
  }

  // Conditional blocks setting specific references, dtref, class, verified based on problem size
  // These assignments are sequential and depend on global problem configuration (grid_points, no_time_steps).
  if (grid_points[0] == 12 && grid_points[1] == 12 && grid_points[2] == 12 && no_time_steps == 60) {
    *class = 'S'; dtref = 1.0e-2;
    xcrref[0] = 1.7034283709541311e-01; xcrref[1] = 1.2975252070034097e-02; xcrref[2] = 3.2527926989486055e-02; xcrref[3] = 2.6436421275166801e-02; xcrref[4] = 1.9211784131744430e-01;
    xceref[0] = 4.9976913345811579e-04; xceref[1] = 4.5195666782961927e-05; xceref[2] = 7.3973765172921357e-05; xceref[3] = 7.3821238632439731e-05; xceref[4] = 8.9269630987491446e-04;
  } else if (grid_points[0] == 24 && grid_points[1] == 24 && grid_points[2] == 24 && no_time_steps == 200) {
    *class = 'W'; dtref = 0.8e-3;
    xcrref[0] = 0.1125590409344e+03; xcrref[1] = 0.1180007595731e+02; xcrref[2] = 0.2710329767846e+02; xcrref[3] = 0.2469174937669e+02; xcrref[4] = 0.2638427874317e+03;
    xceref[0] = 0.4419655736008e+01; xceref[1] = 0.4638531260002e+00; xceref[2] = 1.011551749967e+00; xceref[3] = 0.9235878729944e+00; xceref[4] = 1.018045837718e+01;
  } else if (grid_points[0] == 64 && grid_points[1] == 64 && grid_points[2] == 64 && no_time_steps == 200) {
    *class = 'A'; dtref = 0.8e-3;
    xcrref[0] = 1.0806346714637264e+02; xcrref[1] = 1.1319730901220813e+01; xcrref[2] = 2.5974354511582465e+01; xcrref[3] = 2.3665622544678910e+01; xcrref[4] = 2.5278963211748344e+02;
    xceref[0] = 4.2348416040525025e+00; xceref[1] = 4.4390282496995698e-01; xceref[2] = 9.6692480136345650e-01; xceref[3] = 8.8302063039765474e-01; xceref[4] = 9.7379901770829278e+00;
  } else if (grid_points[0] == 102 && grid_points[1] == 102 && grid_points[2] == 102 && no_time_steps == 200) {
    *class = 'B'; dtref = 3.0e-4;
    xcrref[0] = 1.4233597229287254e+03; xcrref[1] = 9.9330522590150238e+01; xcrref[2] = 3.5646025644535285e+02; xcrref[3] = 3.2485447959084092e+02; xcrref[4] = 3.2707541254659363e+03;
    xceref[0] = 5.2969847140936856e+01; xceref[1] = 4.4632896115670668e+00; xceref[2] = 1.3122573342210174e+01; xceref[3] = 1.2006925323559144e+01; xceref[4] = 1.2459576151035986e+02;
  } else if (grid_points[0] == 162 && grid_points[1] == 162 && grid_points[2] == 162 && no_time_steps == 200) {
    *class = 'C'; dtref = 1.0e-4;
    xcrref[0] = 0.62398116551764615e+04; xcrref[1] = 0.50793239190423964e+03; xcrref[2] = 0.15423530093013596e+04; xcrref[3] = 0.13302387929291190e+04; xcrref[4] = 0.11604087428436455e+05;
    xceref[0] = 0.16462008369091265e+03; xceref[1] = 0.11497107903824313e+02; xceref[2] = 0.41207446207461508e+02; xceref[3] = 0.37087651059694167e+02; xceref[4] = 0.36211053051841265e+03;
  } else {
    // If class remains 'U' after checking known configurations,
    // set verified to FALSE as no reference values are available for comparison.
    *verified = FALSE;
  }

  // Phase 3: Calculate differences (Potentially Parallel)
  // This loop calculates element-wise differences and stores them locally.
  // It depends on results from Phase 1 (xcr, xce) and Phase 2 (xcrref, xceref).
  // Iterations are independent.
  for (m = 0; m < 5; m++) {
    xcrdif[m] = fabs((xcr[m] - xcrref[m]) / xcrref[m]);
    xcedif[m] = fabs((xce[m] - xceref[m]) / xceref[m]);
  }

  // Phase 4: Determine overall verification status (Sequential)
  // This phase checks against epsilon and dtref, updating the shared *verified status.
  // It must be sequential due to the conditional writes to *verified and *class,
  // and the order of checks (DT mismatch causes class='U', which affects subsequent checks).

  // Check DT mismatch first. If class is known, compare dt.
  if (*class != 'U') {
    printf(" Verification being performed for class %1c\n", *class);
    printf(" accuracy setting for epsilon = %20.13e\n", epsilon);
    if (fabs(dt - dtref) > epsilon) {
      *verified = FALSE;
      // Original code sets class to 'U' upon DT mismatch, indicating verification cannot proceed
      // against references, effectively failing verification.
      *class = 'U';
      printf(" DT does not match the reference value of %15.8e\n", dtref);
    }
  } else {
     // If class was 'U' from Phase 2 (unknown problem size), print this.
     // *verified is already FALSE from Phase 2 in this case.
     printf(" Unknown class\n");
  }

  // Now check the difference norms against epsilon.
  // This only makes sense if a known class was successfully determined (*class != 'U').
  // We update *verified if any difference exceeds epsilon.
  // Note: *verified might already be FALSE from Phase 2 or the DT check.
  // Subsequent writes of FALSE don't change the overall failed status.
  if (*class != 'U') {
      // Check residual norm differences (xcrdif)
      for (m = 0; m < 5; m++) {
          if (xcrdif[m] > epsilon) {
              *verified = FALSE; // Found a failure
              break; // No need to check further xcrdif values
          }
      }
      // If not already failed, check solution error norm differences (xcedif)
      if (*verified == TRUE) {
          for (m = 0; m < 5; m++) {
              if (xcedif[m] > epsilon) {
                  *verified = FALSE; // Found a failure
                  break; // No need to check further xcedif values
              }
          }
      }
  }
  // At this point, *verified holds the final verification status based on all checks.

  // Phase 5: Printing results (Sequential)
  // This phase performs I/O and prints the results based on the calculated differences
  // and the final verification status determined in Phase 4.

  // Print header for residual norm comparison
  if (*class != 'U') {
    printf(" Comparison of RMS-norms of residual\n");
  } else {
    printf(" RMS-norms of residual\n");
  }

  // Print residual norms and comparison results - Loop 6
  // Reads local xcr, xcrref, xcrdif and shared *class. Sequential printing.
  for (m = 0; m < 5; m++) {
    if (*class == 'U') {
      // If class is 'U' (unknown or failed DT check), just print computed norm
      printf("          %2d%20.13e\n", m, xcr[m]);
    } else {
      // If class is known, print comparison results
      if (xcrdif[m] > epsilon) {
        printf(" FAILURE: %2d%20.13e%20.13e%20.13e\n",
               m, xcr[m], xcrref[m], xcrdif[m]);
      } else {
        printf("          %2d%20.13e%20.13e%20.13e\n",
               m, xcr[m], xcrref[m], xcrdif[m]);
      }
    }
  }

  // Print header for solution error norm comparison
  if (*class != 'U') {
    printf(" Comparison of RMS-norms of solution error\n");
  } else {
    printf(" RMS-norms of solution error\n");
  }

  // Print solution error norms and comparison results - Loop 7
  // Reads local xce, xceref, xcedif and shared *class. Sequential printing.
  for (m = 0; m < 5; m++) {
    if (*class == 'U') {
      // If class is 'U', just print computed norm
      printf("          %2d%20.13e\n", m, xce[m]);
    } else {
      // If class is known, print comparison results
      if (xcedif[m] > epsilon) {
        printf(" FAILURE: %2d%20.13e%20.13e%20.13e\n",
               m, xce[m], xceref[m], xcedif[m]);
      } else {
        printf("          %2d%20.13e%20.13e%20.13e\n",
               m, xce[m], xceref[m], xcedif[m]);
      }
    }
  }

  // Final verification status print based on the final *verified value
  if (*class == 'U') { // Class 'U' implies verification couldn't be performed or failed DT check
    printf(" No reference values provided\n");
    printf(" No verification performed\n");
  } else if (*verified == TRUE) { // Class is known and all checks passed
    printf(" Verification Successful\n");
  } else { // Class is known, but at least one difference check failed
    printf(" Verification failed\n");
  }
}