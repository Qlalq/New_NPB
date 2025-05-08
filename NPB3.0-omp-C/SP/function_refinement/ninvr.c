#include <stdio.h> // Assuming standard I/O might be needed elsewhere or for debugging
// Add other necessary includes based on original context (e.g., math.h for 0.5, although 0.5 is literal)

// Assuming global or externally defined variables:
// double bt;
// double rhs[5][grid_points[0]][grid_points[1]][grid_points[2]]; // Conceptual declaration
// int grid_points[3]; // Assuming grid dimensions are stored here

static void ninvr(void) {
  int i, j, k;

  // The original structure already provides independent iterations,
  // meaning calculations for rhs[...][i][j][k] only depend on
  // values from rhs[...][i][j][k] itself *before* modification within the same iteration.
  // This is already well-structured for parallelization without complex data dependencies.
  // To make this explicit and potentially aid compilers, we can structure the inner loop
  // to clearly separate reads, calculations, and writes using local temporaries per cell.

  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {

        // 1. Read input values for the current cell (i, j, k)
        double r_in[5];
        r_in[0] = rhs[0][i][j][k];
        r_in[1] = rhs[1][i][j][k];
        r_in[2] = rhs[2][i][j][k];
        r_in[3] = rhs[3][i][j][k];
        r_in[4] = rhs[4][i][j][k];

        // 2. Perform calculations using only the input values and constants (bt)
        double t1 = bt * r_in[2];
        double t2 = 0.5 * ( r_in[3] + r_in[4] );

        // 3. Compute output values for the current cell (i, j, k)
        double r_out[5];
        r_out[0] = -r_in[1];
        r_out[1] = r_in[0];
        r_out[2] = bt * ( r_in[3] - r_in[4] );
        r_out[3] = -t1 + t2;
        r_out[4] = t1 + t2;

        // 4. Write output values back to the array for the current cell (i, j, k)
        // These writes are to unique memory locations for each (i, j, k) iteration,
        // avoiding write-after-read or write-after-write conflicts between iterations.
        rhs[0][i][j][k] = r_out[0];
        rhs[1][i][j][k] = r_out[1];
        rhs[2][i][j][k] = r_out[2];
        rhs[3][i][j][k] = r_out[3];
        rhs[4][i][j][k] = r_out[4];
      }
    }
  }
}