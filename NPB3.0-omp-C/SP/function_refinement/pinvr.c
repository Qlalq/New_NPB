#include <stdio.h> // Assuming standard I/O might be needed elsewhere in a real application
#include <stdlib.h> // Needed for potential dynamic allocation if VLAs are not preferred or supported
// Include other necessary headers for rhs, bt, grid_points if they are part of a larger structure

// Assume the following global or externally defined variables are accessible:
// double rhs[5][GRID_DIM0][GRID_DIM1][GRID_DIM2]; // Replace GRID_DIM with actual sizes or use pointers
// double bt;
// int grid_points[3]; // grid_points[0], grid_points[1], grid_points[2] define the dimensions

// Let's assume the dimensions are available, perhaps defined by constants or accessible variables.
// For clarity in the temporary array declaration, let's assume grid_points is available and large enough.
// We will use Variable Length Arrays (VLAs) for the temporary storage, which is a C99 feature.
// If C99 is not available, dynamic allocation would be necessary.

static void pinvr(void) {
  int i, j, k;
  double r1, r2, r3, r4, r5, t1, t2;

  // Ensure dimensions are valid for the loop bounds (at least 3 in each direction)
  // In a real application, add checks:
  // if (grid_points[0] < 3 || grid_points[1] < 3 || grid_points[2] < 3) { /* handle error */ }

  // Temporary array to store the computed values before writing back to rhs.
  // This decouples the computation (reading original values, writing new values)
  // from the final update step, removing any read-after-write or write-after-read
  // dependencies within a single (i, j, k) block during the computation phase
  // if they were present (though in this specific case, the original code
  // already reads all needed values before writing).
  // This structure makes the two distinct phases (compute, update) explicit
  // and makes both loops trivially parallelizable.

  // Note: Using VLA requires C99 or C11. For C89 or dynamic dimensions at runtime
  // without VLA, dynamic allocation using malloc would be required.
  // Example using VLA:
  double temp_rhs[5][grid_points[0]][grid_points[1]][grid_points[2]];

  // --- Pass 1: Compute new values and store in temporary array ---
  // This loop can be parallelized directly using OpenMP's #pragma omp for
  // (or #pragma omp parallel for). No loop-carried dependencies or write conflicts
  // on temp_rhs for different (i, j, k) iterations.
  // Workload per iteration is constant, ensuring good load balance.
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	// Reads from original rhs (thread-safe as reads are independent per iteration)
	r1 = rhs[0][i][j][k];
	r2 = rhs[1][i][j][k];
	r3 = rhs[2][i][j][k];
	r4 = rhs[3][i][j][k];
	r5 = rhs[4][i][j][k];

	// Computations (thread-local variables)
	t1 = bt * r1;
	t2 = 0.5 * ( r4 + r5 );

	// Writes to temporary storage (thread-safe as each iteration writes to unique temp_rhs location)
	temp_rhs[0][i][j][k] =  bt * ( r4 - r5 );
	temp_rhs[1][i][j][k] = -r3;
	temp_rhs[2][i][j][k] =  r2;
	temp_rhs[3][i][j][k] = -t1 + t2;
	temp_rhs[4][i][j][k] =  t1 + t2;
      }
    }
  }

  // --- Pass 2: Copy computed values from temporary array back to original rhs ---
  // This loop can also be parallelized directly using OpenMP's #pragma omp for.
  // Reads from temp_rhs and writes to rhs are independent for different (i, j, k) iterations.
  // Workload per iteration is constant (simple copy), ensuring good load balance.
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	// Writes to original rhs (thread-safe as each iteration writes to unique rhs location)
	rhs[0][i][j][k] = temp_rhs[0][i][j][k];
	rhs[1][i][j][k] = temp_rhs[1][i][j][k];
	rhs[2][i][j][k] = temp_rhs[2][i][j][k];
	rhs[3][i][j][k] = temp_rhs[3][i][j][k];
	rhs[4][i][j][k] = temp_rhs[4][i][j][k];
      }
    }
  }

  // If dynamic allocation was used for temp_rhs, it would need to be freed here.
}