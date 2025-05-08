#include <stdio.h>

// Assuming grid_points is a global or accessible array, e.g.:
// extern int grid_points[3];

// Assuming definitions for BB, CC, AA constants, e.g.:
// #define BB 0
// #define CC 1
// #define AA 2

// Assuming declarations for lhs and rhs arrays, e.g.:
// extern double lhs[GRID_SIZE_X][GRID_SIZE_Y][GRID_SIZE_Z][NUM_MATRIX_ELEMENTS];
// extern double rhs[GRID_SIZE_X][GRID_SIZE_Y][GRID_SIZE_Z];

// Assuming definitions for helper functions:
// void binvcrhs(double* b, double* c, double* r);
// void matvec_sub(double* a, double* r_in, double* r_out);
// void matmul_sub(double* a, double* c_in, double* b_out);
// void binvrhs(double* b, double* r);


static void y_solve_cell(void) {
  int i, j, k;
  int jsize;
  int imin, imax, kmin, kmax;

  // jsize is the upper bound for j in the main loop iteration
  // and also the index for the final boundary plane.
  // The j loop iterates from 1 up to jsize-1.
  jsize = grid_points[1]-1;

  // Define constant loop bounds for dimensions i and k.
  // These loops iterate over the interior points (excluding boundaries at 0 and size-1).
  imin = 1;
  imax = grid_points[0]-1; // Loop goes up to imax-1
  kmin = 1;
  kmax = grid_points[2]-1; // Loop goes up to kmax-1


  // Part 1: Processing the boundary plane j = 0.
  // The computations for different (i, k) points are independent of each other.
  // This nested loop block is a candidate for OpenMP parallelization over i and k.
  for (i = imin; i < imax; i++) {
    for (k = kmin; k < kmax; k++) {
      // binvcrhs updates lhs[i][0][k][CC] and rhs[i][0][k]
      binvcrhs( lhs[i][0][k][BB],
		lhs[i][0][k][CC],
		rhs[i][0][k] );
    }
  }

  // Part 2: Marching steps along the j dimension (from j=1 to jsize-1).
  // The outer loop over 'j' must be sequential due to loop-carried dependencies:
  // the computation at step 'j' depends on results from step 'j-1'.
  // However, for a *fixed* value of 'j', the computations for different (i, k) points
  // are independent of each other.
  // The nested loops over i and k are candidates for OpenMP parallelization inside the j loop.
  for (j = 1; j < jsize; j++) { // j goes from 1 up to (grid_points[1]-1)-1 = grid_points[1]-2
    for (i = imin; i < imax; i++) {
      for (k = kmin; k < kmax; k++) {
	// matvec_sub uses rhs[i][j-1][k], writes rhs[i][j][k]
	matvec_sub(lhs[i][j][k][AA],
		   rhs[i][j-1][k], rhs[i][j][k]);
	// matmul_sub uses lhs[i][j-1][k][CC], writes lhs[i][j][k][BB]
	matmul_sub(lhs[i][j][k][AA],
		   lhs[i][j-1][k][CC],
		   lhs[i][j][k][BB]);
	// binvcrhs uses/writes lhs[i][j][k][CC] and rhs[i][j][k]
	binvcrhs( lhs[i][j][k][BB],
		  lhs[i][j][k][CC],
		  rhs[i][j][k] );
      }
    }
  }

  // Part 3: Processing the boundary plane j = jsize (which is grid_points[1]-1).
  // Similar to Part 1, computations for different (i, k) points are independent.
  // This nested loop block is a candidate for OpenMP parallelization over i and k.
  // This step uses results from the last j iteration (jsize-1).
  for (i = imin; i < imax; i++) {
    for (k = kmin; k < kmax; k++) {
      // matvec_sub uses rhs[i][jsize-1][k], writes rhs[i][jsize][k]
      matvec_sub(lhs[i][jsize][k][AA],
		 rhs[i][jsize-1][k], rhs[i][jsize][k]);
      // matmul_sub uses lhs[i][jsize-1][k][CC], writes lhs[i][jsize][k][BB]
      matmul_sub(lhs[i][jsize][k][AA],
		 lhs[i][jsize-1][k][CC],
		 lhs[i][jsize][k][BB]);
      // binvrhs uses lhs[i][jsize][k][BB], writes rhs[i][jsize][k]
      binvrhs( lhs[i][jsize][k][BB],
	       rhs[i][jsize][k] );
    }
  }
}