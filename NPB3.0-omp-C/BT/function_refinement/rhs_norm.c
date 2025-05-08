#include <math.h> // For sqrt

// Assume grid_points and rhs are globally accessible or defined elsewhere
// based on the function signature and usage in the original code.
// Example declaration assuming max sizes or external definition:
// extern int grid_points[3];
// extern double rhs[MAX_DIM1][MAX_DIM2][MAX_DIM3][5];


static void rhs_norm_refined(double rms[5]) {
  int i, j, k, d;

  // Initialize rms. This loop is independent across m=0..4.
  // It can be parallelized or run sequentially as it's small.
  for (int m = 0; m < 5; m++) {
    rms[m] = 0.0;
  }

  // Refactored structure to explicitly separate computation for each 'm' component (0 to 4).
  // Each block below accumulates into a distinct element of the 'rms' array,
  // making the blocks independent of each other with respect to 'rms' updates.
  // This structure is highly amenable to OpenMP 'sections' for parallelization of the 5 blocks.

  // Accumulation for rms[0] (m=0)
  int m = 0; // Fix m for this independent block
  for (i = 1; i < grid_points[0]-1; i++) {
    // The nested i, j, k loops process a large 3D grid.
    // These inner loops represent a significant workload and can be parallelized
    // within this block (e.g., using OpenMP parallel for with reduction on rms[m]).
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
        double add = rhs[i][j][k][m];
        rms[m] = rms[m] + add*add; // Accumulating into rms[0]
      }
    }
  }

  // Accumulation for rms[1] (m=1)
  m = 1; // Fix m for this independent block
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
        double add = rhs[i][j][k][m];
        rms[m] = rms[m] + add*add; // Accumulating into rms[1]
      }
    }
  }

  // Accumulation for rms[2] (m=2)
  m = 2; // Fix m for this independent block
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
        double add = rhs[i][j][k][m];
        rms[m] = rms[m] + add*add; // Accumulating into rms[2]
      }
    }
  }

  // Accumulation for rms[3] (m=3)
  m = 3; // Fix m for this independent block
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
        double add = rhs[i][j][k][m];
        rms[m] = rms[m] + add*add; // Accumulating into rms[3]
      }
    }
  }

  // Accumulation for rms[4] (m=4)
  m = 4; // Fix m for this independent block
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
        double add = rhs[i][j][k][m];
        rms[m] = rms[m] + add*add; // Accumulating into rms[4]
      }
    }
  }

  // Finalization loop. Independent across m=0..4.
  // This loop can be parallelized after all accumulations are complete.
  for (m = 0; m < 5; m++) {
    // Inner loop over d is small (3 iterations)
    for (d = 0; d <= 2; d++) {
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    }
    rms[m] = sqrt(rms[m]);
  }
}