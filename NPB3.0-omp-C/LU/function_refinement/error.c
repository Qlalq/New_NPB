#include <math.h> // Required for sqrt

// Assuming errnm is a global or file-static array accessible here, e.g.:
// double errnm[5];
// Assuming ist, iend, jst, jend, nz, nx0, ny0, nz0 are accessible constants or variables.
// Assuming u is a global or file-static multi-dimensional array accessible here, e.g.:
// double u[...][...][...][5];
// Assuming exact is a function with signature void exact(int, int, int, double*);

static void error_refined(void) {
  int i, j, k, m;
  double tmp;
  double u000ijk[5];

  // --- Step 1: Initialization ---
  // Initialize the final result accumulator.
  // This loop is sequential and typically fast.
  for (m = 0; m < 5; m++) {
    errnm[m] = 0.0;
  }

  // --- Step 2: Main computation with decoupled accumulation ---
  // Identify the range of the primary parallelizable loop ('i' in this case).
  // If iend < ist, the loop range is empty. Handle this case to avoid VLA size issues.
  int i_range_size = 0;
  if (iend >= ist) {
      i_range_size = iend - ist + 1;
  }

  // Create temporary storage to hold partial sums for each 'i' iteration.
  // This decouples the accumulation dependency between iterations of the 'i' loop.
  // Using Variable Length Array (VLA) - requires C99 or later.
  // If C89 is required, this would need dynamic allocation (malloc).
  double partial_sums[i_range_size][5];

  // Initialize the temporary storage. This loop can potentially be parallelized as well.
  for (int ii = 0; ii < i_range_size; ++ii) {
      for (m = 0; m < 5; ++m) {
          partial_sums[ii][m] = 0.0;
      }
  }

  // Compute the partial sums for each 'i' iteration.
  // This loop nest is the primary target for OpenMP parallelization over 'i'.
  // Iterations of this 'i' loop write to distinct locations in 'partial_sums',
  // removing the loop-carried dependency that was on 'errnm' in the original code.
  for (i = ist; i <= iend; i++) {
    // Map 'i' index (ist..iend) to 0-based index for partial_sums (0..i_range_size-1).
    int i_idx = i - ist;

    // The inner loops iterate over j and k for the current i.
    for (j = jst; j <= jend; j++) {
      for (k = 1; k <= nz - 2; k++) {
        // Call exact function - assumes it reads i, j, k and writes to u000ijk.
        // u000ijk is local to this (i,j,k) block, or thread-local in parallel context.
        exact(i, j, k, u000ijk); // iglob and jglob are just i and j, removed them for simplicity

        // Accumulate sum of squares into the temporary array, unique to this 'i'.
        for (m = 0; m < 5; m++) {
          tmp = (u000ijk[m] - u[i][j][k][m]); // u is assumed read-only
          partial_sums[i_idx][m] += tmp * tmp; // Accumulate into temporary array
        }
      }
    }
  }
  // After this loop nest, partial_sums[ii][m] holds the sum of squares for m,
  // accumulated over j and k, for the specific i = ist + ii.

  // --- Step 3: Combine partial sums (Reduction step) ---
  // Sum the contributions from all 'i' iterations (stored in partial_sums) into the final 'errnm'.
  // This is an explicit reduction step.
  // This loop can be parallelized using OpenMP's reduction clause on errnm[m]
  // (e.g., over the outer 'm' loop or inner 'ii' loop).
  for (m = 0; m < 5; m++) {
    for (int ii = 0; ii < i_range_size; ++ii) {
      errnm[m] += partial_sums[ii][m]; // Accumulate from temporary storage
    }
  }
  // Now errnm[m] contains the total sum of squares over all i, j, k.

  // --- Step 4: Final calculation ---
  // Perform the final calculation on the accumulated sum.
  // This loop depends on the completed values of errnm and runs sequentially after the reduction.
  // Ensure the denominator is not zero to prevent division by zero or NaN.
  double denominator = (double)(nx0 - 2) * (ny0 - 2) * (nz0 - 2);
  if (denominator == 0.0) {
      // Handle error case or set errnm to infinity/NaN as appropriate
      // For this example, we might just let it result in Inf/NaN if input sizes are invalid
  }

  for (m = 0; m < 5; m++) {
    // errnm[m] should be non-negative as it's a sum of squares.
    if (errnm[m] < 0.0) errnm[m] = 0.0; // Guard against potential negative due to precision issues
    errnm[m] = sqrt(errnm[m] / denominator);
  }
}