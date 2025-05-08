#include <stdio.h>
#include <math.h> // Required for sqrt

// Assuming grid_points and rhs are defined globally or accessible in scope
// For compilation purposes, let's assume minimal definitions if they are not provided globally.
// You should replace these with the actual definitions from your project.
// extern int grid_points[3]; // Example: extern int grid_points[3] = {10, 10, 10};
// extern double *****rhs;     // Example: This would be a complex declaration depending on the actual dimensions and layout.
                             // Assuming rhs is double rhs[5][grid_points[0]-1][grid_points[1]-1][grid_points[2]-1]
                             // Accessing rhs[m][i][j][k] implies dimensions [5][dim0][dim1][dim2]

static void rhs_norm_refined(double rms[5]) {
  int i, j, k, d, m;
  double add;

  // Allocate a temporary array to store the raw sum of squares for each component.
  // This decouples the pure accumulation phase from the final post-processing steps (division, sqrt).
  // By accumulating into a temporary array first, the dependency on the final 'rms' array
  // during the main computation loop is removed. The 'sum_of_squares' array becomes
  // the explicit target for reduction operations when parallelizing.
  double sum_of_squares[5];

  // Initialize the temporary accumulation array.
  // This loop is simple and fast. It can be sequential or parallelized if needed.
  for (m = 0; m < 5; m++) {
    sum_of_squares[m] = 0.0;
  }

  // Main accumulation loop nest.
  // This is the computationally intensive part.
  // The loops over i, j, and k iterate through the grid points (excluding boundaries used in the original).
  // Iterations across the (i, j, k) space read from 'rhs' independently for a fixed 'm'.
  // The accumulation into 'sum_of_squares[m]' is a reduction operation across these (i, j, k) iterations.
  // Structuring with i, j, k as outer loops makes them the primary candidates for
  // parallel decomposition (e.g., using OpenMP parallel for with collapse).
  // The inner loop over m iterates through the 5 components and performs the accumulation
  // for the current grid point (i, j, k). Keeping this loop inner might be beneficial
  // for cache locality when accessing rhs[m][i][j][k].
  for (i = 0; i <= grid_points[0]-2; i++) {
    for (j = 0; j <= grid_points[1]-2; j++) {
      for (k = 0; k <= grid_points[2]-2; k++) {
	    for (m = 0; m < 5; m++) {
	      add = rhs[m][i][j][k];
	      // Accumulate the square of the value into the corresponding component's sum
	      // in the temporary array. This is the operation that requires reduction
	      // handling if the i, j, k loops are parallelized.
	      sum_of_squares[m] = sum_of_squares[m] + add*add;
	    }
      }
    }
  }

  // Final calculation steps: calculate the mean square and take the square root.
  // This loop iterates through the 5 components.
  // The calculation for each component 'rms[m]' is independent of other components 'rms[m']'.
  // It uses the pre-calculated 'sum_of_squares[m]'.
  // This loop is clearly parallelizable over 'm' as there are no dependencies between iterations.
  for (m = 0; m < 5; m++) {
    // Start with the total sum of squares accumulated for component m.
    rms[m] = sum_of_squares[m];

    // Divide by the number of points in each dimension (grid_points[d]-2) to get the mean square.
    // This loop is short (3 iterations) and calculates the average over the computational domain size.
    for (d = 0; d < 3; d++) {
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    }
    // Take the square root to get the root mean square for component m.
    rms[m] = sqrt(rms[m]);
  }
}