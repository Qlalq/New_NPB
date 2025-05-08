#include <math.h> // Required for sqrt
#include <stdio.h> // Assuming for grid_points or potential debugging, can be removed if not used

// Assume grid_points, dnxm1, dnym1, dnzm1, u, and exact_solution are defined elsewhere
// Example external declarations (not part of the function itself):
// extern int grid_points[3];
// extern double dnxm1, dnym1, dnzm1;
// extern double u[5][ANY_SIZE_X][ANY_SIZE_Y][ANY_SIZE_Z]; // Replace ANY_SIZE with actual dimensions
// extern void exact_solution(double xi, double eta, double zeta, double u_exact[5]);


static void error_norm(double rms[5]) {
  int i, j, k, m; // d variable scope is moved

  // Variable declarations for use within the main loops
  double xi, eta, zeta, u_exact[5], add;

  // 1. Initialization:
  // This loop initializes the reduction variables. Can be parallelized
  // if rms is thread-private before the reduction, but simple sequential
  // initialization is often fine if the main loop is parallel.
  for (m = 0; m < 5; m++) {
    rms[m] = 0.0;
  }

  // 2. Calculation of sum of squares (Main Reduction Loop):
  // This triple nested loop over i, j, k is the primary target for parallelization.
  // The inner loop accumulates sum of squares into rms[m].
  // This is a reduction operation on the rms array.
  // The current structure is a correct sequential reduction.
  // For OpenMP, this loop nest (or the outermost loop of it) would be parallelized
  // with a reduction clause on `rms`. No structural change is needed in the
  // sequential logic to make it *represent* a reduction; it already does.
  // Variables xi, eta, zeta, u_exact, and add are loop-local and suitable for privatization.
  for (i = 0; i <= grid_points[0]-1; i++) {
    xi = (double)i * dnxm1;
    for (j = 0; j <= grid_points[1]-1; j++) {
      eta = (double)j * dnym1;
      for (k = 0; k <= grid_points[2]-1; k++) {
	zeta = (double)k * dnzm1;

	// exact_solution depends on current (xi, eta, zeta)
	exact_solution(xi, eta, zeta, u_exact);

	// Inner loop accumulates squared differences for each of the 5 components (m).
	// This is the core reduction update: rms[m] += add*add;
	for (m = 0; m < 5; m++) {
	  add = u[m][i][j][k] - u_exact[m];
	  rms[m] = rms[m] + add*add; // Accumulation (Reduction)
	}
      }
    }
  }

  // 3. Final calculation (division and sqrt):
  // The original code had a loop-carried dependency for rms[m] in the d loop
  // (rms[m] = rms[m] / (grid_points[d]-2)).
  // Decouple this by calculating the total divisor outside the loop.
  double total_divisor = 1.0;
  int d; // d is now local to this calculation
  for (d = 0; d < 3; d++) {
     // Assumes grid_points[d]-2 is the factor for each dimension
     total_divisor *= (double)(grid_points[d]-2);
  }

  // Now, the loop over m is completely independent, operating on the final
  // accumulated rms values. This loop is trivially parallelizable.
  for (m = 0; m < 5; m++) {
    // Avoid division by zero, though for valid grid sizes > 2 this shouldn't happen.
    if (total_divisor != 0.0) {
       rms[m] = rms[m] / total_divisor;
    } else {
       // Handle unexpected case: set to 0 or report error
       rms[m] = 0.0;
    }
    // Final square root calculation
    rms[m] = sqrt(rms[m]);
  }
}