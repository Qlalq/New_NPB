#include <math.h> // Required for sqrt

static void error_norm(double rms[5]) {
  int i, j, k, m, d;
  double xi, eta, zeta, u_exact[5], add;
  double sum_m; // Temporary accumulator for reduction for each m

  // Initialization of rms is implicitly handled by initializing sum_m
  // within the loop over m below. If an explicit zeroing of the
  // final 'rms' array before the loop is desired, it could be added here,
  // but the current structure overwrites it anyway.

  // Main Computation: Restructure to make accumulation for each rms[m] independent.
  // Loop over m (0 to 4) - this loop's iterations are independent, as each iteration
  // calculates the sum for a specific rms[m] element without depending on other rms elements.
  // This structure allows for potential parallelization of this outer loop over m.
  for (m = 0; m < 5; m++) {
    // Initialize reduction variable for this specific rms[m] element.
    // This variable is local to the conceptual work unit for m.
    sum_m = 0.0;

    // Loop over the 3D grid points (i, j, k).
    // The sum calculation within these loops for a fixed m is a reduction.
    // These loops represent the bulk of the work and can be parallelized
    // using a reduction clause on sum_m.
    for (i = 0; i < grid_points[0]; i++) {
      xi = (double)i * dnxm1;
      for (j = 0; j < grid_points[1]; j++) {
        eta = (double)j * dnym1;
        for (k = 0; k < grid_points[2]; k++) {
          zeta = (double)k * dnzm1;

          // Calculate exact solution at (i,j,k) - result is u_exact[0..4].
          // This call depends on i,j,k but is independent of m at this point.
          exact_solution(xi, eta, zeta, u_exact);

          // Calculate squared difference for the current 'm'.
          add = u[i][j][k][m] - u_exact[m];

          // Accumulate the squared difference for the current 'm' into the sum_m.
          // This is the core reduction operation that happens across the i,j,k space.
          sum_m += add * add;
        }
      }
    }
    // After iterating through all i,j,k points for a specific m, store the accumulated sum
    // into the corresponding element of the rms array. This write is safe because each
    // m iteration (or thread handling an m iteration) writes to a distinct rms[m].
    rms[m] = sum_m;
  }

  // Final Calculation: Normalize and take square root.
  // Loop over m (0 to 4) - this loop's iterations are independent.
  // This loop can be parallelized over m.
  for (m = 0; m < 5; m++) {
    // Normalize by dividing by (grid_points[d]-2) for d=0,1,2.
    // This inner loop has a loop-carried dependency on rms[m], where rms[m] is updated
    // based on its previous value within the d loop. However, this loop is small (3 iterations).
    for (d = 0; d <= 2; d++) {
      rms[m] = rms[m] / (double)(grid_points[d] - 2);
    }
    // Take the square root of the normalized sum. Depends on the value of rms[m] after the d loop.
    rms[m] = sqrt(rms[m]);
  }
}