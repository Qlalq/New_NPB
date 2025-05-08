#include <stdio.h> // Assuming stdio.h might be needed based on original code's context (e.g., for printing, though not shown in the function)

// Assumed global variables based on the original code's usage.
// These would need proper declarations in a header file or before the function.
// Example declarations (actual types/dimensions depend on the full program):
// extern int grid_points[3];
// extern double lhs[5][XDIM][YDIM][ZDIM]; // Replace XDIM, YDIM, ZDIM with actual sizes derived from grid_points
// extern double rhs[5][XDIM][YDIM][ZDIM]; // Replace XDIM, YDIM, ZDIM with actual sizes derived from grid_points

// Assuming lhsx and ninvr are defined elsewhere
extern void lhsx(void);
extern void ninvr(void);

// Refined function for improved parallelization potential
// The core loops over 'i' have dependencies that prevent simple parallelization.
// The loops over 'j', 'k', and 'm' are independent for a fixed 'i' and are
// the primary targets for OpenMP parallelization (e.g., using #pragma omp parallel for).
static void x_solve_refined(void) {
  int i, j, k, m;
  double fac1, fac2;

  lhsx(); // Assume lhsx is a setup phase, likely sequential

  //---------------------------------------------------------------------
  // Forward elimination
  // Loop over 'i' has dependencies (i depends on i-1, i-2 in the backward solve
  // part, and updates i+1, i+2 based on i in this forward sweep).
  // The loops over 'j' and 'k' are independent for a fixed 'i'.
  // The inner 'm' loop for rhs updates is also independent.
  // Parallelization target: Loops over j, k, and/or m.
  //---------------------------------------------------------------------

  // Block 1 & 3: i = 0 to grid_points[0]-3
  // Split based on 'm' range as in original code due to different 'n' offset
  int i_max1 = grid_points[0] - 3;

  // m = 0 to 2 section
  int n_0 = 0; // n is fixed at 0 for m=0..2
  for (i = 0; i <= i_max1; ++i) { // Sequential loop due to i-dependency
    int i1 = i + 1;
    int i2 = i + 2;

    // These loops over j, k are independent and can be parallelized.
    // For example: #pragma omp parallel for collapse(2) private(j, k, fac1, m)
    for (j = 1; j <= grid_points[1]-2; ++j) {
      for (k = 1; k <= grid_points[2]-2; ++k) {
	fac1 = 1./lhs[n_0+2][i][j][k];
	lhs[n_0+3][i][j][k]   = fac1*lhs[n_0+3][i][j][k];
	lhs[n_0+4][i][j][k]   = fac1*lhs[n_0+4][i][j][k];
	for (m = 0; m < 3; ++m) { // This m loop is also independent
	  rhs[m][i][j][k] = fac1*rhs[m][i][j][k];
	}

	lhs[n_0+2][i1][j][k] = lhs[n_0+2][i1][j][k] -
	  lhs[n_0+1][i1][j][k]*lhs[n_0+3][i][j][k];
	lhs[n_0+3][i1][j][k] = lhs[n_0+3][i1][j][k] -
	  lhs[n_0+1][i1][j][k]*lhs[n_0+4][i][j][k];
	for (m = 0; m < 3; ++m) { // This m loop is also independent
	  rhs[m][i1][j][k] = rhs[m][i1][j][k] -
	    lhs[n_0+1][i1][j][k]*rhs[m][i][j][k];
	}

	lhs[n_0+1][i2][j][k] = lhs[n_0+1][i2][j][k] -
	  lhs[n_0+0][i2][j][k]*lhs[n_0+3][i][j][k];
	lhs[n_0+2][i2][j][k] = lhs[n_0+2][i2][j][k] -
	  lhs[n_0+0][i2][j][k]*lhs[n_0+4][i][j][k];
	for (m = 0; m < 3; ++m) { // This m loop is also independent
	  rhs[m][i2][j][k] = rhs[m][i2][j][k] -
	    lhs[n_0+0][i2][j][k]*rhs[m][i][j][k];
	}
      }
    }
  }

  // m = 3 to 4 section
  // The outer m loop is independent and can be parallelized.
  // The inner j, k loops are also independent.
  // Parallelization target: Outer m loop OR loops over j and k.
  for (m = 3; m < 5; ++m) { // Potential parallel loop over m
    int n_m = (m-3+1)*5;

    for (i = 0; i <= i_max1; ++i) { // Sequential loop due to i-dependency
      int i1 = i + 1;
      int i2 = i + 2;

      // These loops over j, k are independent and can be parallelized.
      // For example: #pragma omp parallel for collapse(2) private(j, k, fac1)
      for (j = 1; j <= grid_points[1]-2; ++j) {
	for (k = 1; k <= grid_points[2]-2; ++k) {
	  fac1               = 1./lhs[n_m+2][i][j][k];
	  lhs[n_m+3][i][j][k]   = fac1*lhs[n_m+3][i][j][k];
	  lhs[n_m+4][i][j][k]   = fac1*lhs[n_m+4][i][j][k];
	  rhs[m][i][j][k] = fac1*rhs[m][i][j][k];

	  lhs[n_m+2][i1][j][k] = lhs[n_m+2][i1][j][k] -
	    lhs[n_m+1][i1][j][k]*lhs[n_m+3][i][j][k];
	  lhs[n_m+3][i1][j][k] = lhs[n_m+3][i1][j][k] -
	    lhs[n_m+1][i1][j][k]*lhs[n_m+4][i][j][k];
	  rhs[m][i1][j][k] = rhs[m][i1][j][k] -
	    lhs[n_m+1][i1][j][k]*rhs[m][i][j][k];

	  lhs[n_m+1][i2][j][k] = lhs[n_m+1][i2][j][k] -
	    lhs[n_m+0][i2][j][k]*lhs[n_m+3][i][j][k];
	  lhs[n_m+2][i2][j][k] = lhs[n_m+2][i2][j][k] -
	    lhs[n_m+0][i2][j][k]*lhs[n_m+4][i][j][k];
	  rhs[m][i2][j][k] = rhs[m][i2][j][k] -
	    lhs[n_m+0][i2][j][k]*rhs[m][i][j][k];
	}
      }
    }
  }


  // Block 2 & 4: i = grid_points[0]-2 (handling the last two planes)
  // The loop over i is fixed here, j and k loops are independent.
  // Split based on 'm' range as in original code.
  int i_last_minus_1 = grid_points[0]-2;
  int i_last         = grid_points[0]-1;

  // m = 0 to 2 section
  n_0 = 0; // n is fixed at 0
  // These loops over j, k are independent and can be parallelized.
  // For example: #pragma omp parallel for collapse(2) private(j, k, fac1, fac2, m)
  for (j = 1; j <= grid_points[1]-2; ++j) {
    for (k = 1; k <= grid_points[2]-2; ++k) {
      fac1               = 1.0/lhs[n_0+2][i_last_minus_1][j][k];
      lhs[n_0+3][i_last_minus_1][j][k]   = fac1*lhs[n_0+3][i_last_minus_1][j][k];
      lhs[n_0+4][i_last_minus_1][j][k]   = fac1*lhs[n_0+4][i_last_minus_1][j][k];
      for (m = 0; m < 3; ++m) { // This m loop is independent
	rhs[m][i_last_minus_1][j][k] = fac1*rhs[m][i_last_minus_1][j][k];
      }

      lhs[n_0+2][i_last][j][k] = lhs[n_0+2][i_last][j][k] -
	lhs[n_0+1][i_last][j][k]*lhs[n_0+3][i_last_minus_1][j][k];
      lhs[n_0+3][i_last][j][k] = lhs[n_0+3][i_last][j][k] -
	lhs[n_0+1][i_last][j][k]*lhs[n_0+4][i_last_minus_1][j][k];
      for (m = 0; m < 3; ++m) { // This m loop is independent
	rhs[m][i_last][j][k] = rhs[m][i_last][j][k] -
	  lhs[n_0+1][i_last][j][k]*rhs[m][i_last_minus_1][j][k];
      }

      fac2               = 1./lhs[n_0+2][i_last][j][k];
      for (m = 0; m < 3; ++m) { // This m loop is independent
	rhs[m][i_last][j][k] = fac2*rhs[m][i_last][j][k];
      }
    }
  }

  // m = 3 to 4 section
  // The outer m loop is independent and can be parallelized.
  // The inner j, k loops are also independent.
  // Parallelization target: Outer m loop OR loops over j and k.
  for (m = 3; m < 5; ++m) { // Potential parallel loop over m
    int n_m = (m-3+1)*5;

    // These loops over j, k are independent and can be parallelized.
    // For example: #pragma omp parallel for collapse(2) private(j, k, fac1, fac2)
    for (j = 1; j <= grid_points[1]-2; ++j) {
      for (k = 1; k <= grid_points[2]-2; ++k) {
	fac1               = 1./lhs[n_m+2][i_last_minus_1][j][k];
	lhs[n_m+3][i_last_minus_1][j][k]   = fac1*lhs[n_m+3][i_last_minus_1][j][k];
	lhs[n_m+4][i_last_minus_1][j][k]   = fac1*lhs[n_m+4][i_last_minus_1][j][k];
	rhs[m][i_last_minus_1][j][k] = fac1*rhs[m][i_last_minus_1][j][k];

	lhs[n_m+2][i_last][j][k] = lhs[n_m+2][i_last][j][k] -
	  lhs[n_m+1][i_last][j][k]*lhs[n_m+3][i_last_minus_1][j][k];
	lhs[n_m+3][i_last][j][k] = lhs[n_m+3][i_last][j][k] -
	  lhs[n_m+1][i_last][j][k]*lhs[n_m+4][i_last_minus_1][j][k];
	rhs[m][i_last][j][k] = rhs[m][i_last][j][k] -
	  lhs[n_m+1][i_last][j][k]*rhs[m][i_last_minus_1][j][k];

	fac2               = 1./lhs[n_m+2][i_last][j][k];
	rhs[m][i_last][j][k] = fac2*rhs[m][i_last][j][k];
      }
    }
  }

  //---------------------------------------------------------------------
  // Backward substitution
  // Loop over 'i' is backward and has dependencies (i depends on i+1, i+2).
  // The loops over 'm', 'j', and 'k' are independent for a fixed 'i'.
  // Parallelization target: Loops over m, j, and/or k.
  //---------------------------------------------------------------------

  // Block 5 & 6: i = grid_points[0]-2 (handling the last two planes)
  // Split based on 'm' range. m, j, k loops are independent.
  // Parallelization target: Loops over m, j, and/or k.
  i_last_minus_1 = grid_points[0]-2;
  i_last         = grid_points[0]-1;

  // m = 0 to 2 section
  n_0 = 0; // n is fixed at 0
  // These loops over m, j, k are independent and can be parallelized.
  // For example: #pragma omp parallel for collapse(3) private(m, j, k)
  for (m = 0; m < 3; ++m) {
    for (j = 1; j <= grid_points[1]-2; ++j) {
      for (k = 1; k <= grid_points[2]-2; ++k) {
	rhs[m][i_last_minus_1][j][k] = rhs[m][i_last_minus_1][j][k] -
	  lhs[n_0+3][i_last_minus_1][j][k]*rhs[m][i_last][j][k];
      }
    }
  }

  // m = 3 to 4 section
  // These loops over m, j, k are independent and can be parallelized.
  // For example: #pragma omp parallel for collapse(3) private(m, j, k, n_m)
  for (m = 3; m < 5; ++m) {
    int n_m = (m-3+1)*5;
    for (j = 1; j <= grid_points[1]-2; ++j) {
      for (k = 1; k <= grid_points[2]-2; ++k) {
	rhs[m][i_last_minus_1][j][k] = rhs[m][i_last_minus_1][j][k] -
	  lhs[n_m+3][i_last_minus_1][j][k]*rhs[m][i_last][j][k];
      }
    }
  }

  // Block 7 & 8: i = grid_points[0]-3 down to 0
  // Split based on 'm' range.
  // Loop over 'i' is sequential due to dependencies.
  // The loops over 'm', 'j', and 'k' are independent for a fixed 'i'.
  // Parallelization target: Loops over m, j, and/or k.

  i_max1 = grid_points[0]-3;

  // m = 0 to 2 section
  n_0 = 0; // n is fixed at 0
  for (i = i_max1; i >= 0; --i) { // Sequential loop due to i-dependency
    int i1 = i + 1;
    int i2 = i + 2;

    // These loops over m, j, k are independent and can be parallelized.
    // For example: #pragma omp parallel for collapse(3) private(m, j, k)
    for (m = 0; m < 3; ++m) {
      for (j = 1; j <= grid_points[1]-2; ++j) {
	for (k = 1; k <= grid_points[2]-2; ++k) {
	  rhs[m][i][j][k] = rhs[m][i][j][k] -
	    lhs[n_0+3][i][j][k]*rhs[m][i1][j][k] -
	    lhs[n_0+4][i][j][k]*rhs[m][i2][j][k];
	}
      }
    }
  }

  // m = 3 to 4 section
  // The outer m loop is independent and can be parallelized.
  // The inner j, k loops are also independent.
  // Parallelization target: Outer m loop OR loops over j and k.
  for (m = 3; m < 5; ++m) { // Potential parallel loop over m
    int n_m = (m-3+1)*5;

    for (i = i_max1; i >= 0; --i) { // Sequential loop due to i-dependency
      int i1 = i + 1;
      int i2 = i + 2;

      // These loops over j, k are independent and can be parallelized.
      // For example: #pragma omp parallel for collapse(2) private(j, k)
      for (j = 1; j <= grid_points[1]-2; ++j) {
	for (k = 1; k <= grid_points[2]-2; ++k) {
	  rhs[m][i][j][k] = rhs[m][i][j][k] -
	    lhs[n_m+3][i][j][k]*rhs[m][i1][j][k] -
	    lhs[n_m+4][i][j][k]*rhs[m][i2][j][k];
	}
      }
    }
  }

  ninvr(); // Assume ninvr is a finalization phase, likely sequential
}