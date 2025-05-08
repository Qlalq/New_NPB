#include <math.h> // Required for sqrt function
// Assuming ISIZ1, ISIZ2, ISIZ3 are defined elsewhere, e.g., in a header file
// #include "global_size_definitions.h"

static void l2norm (int nx0, int ny0, int nz0,
		    int ist, int iend,
		    int jst, int jend,
		    double v[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5],
		    double sum[5]) {

  int i, j, k, m;
  double norm_factor; // Variable to store the normalization factor

  // 1. Data Dependencies & Synchronization Overhead:
  // The initialization of the sum array is independent and sequential.
  // No complex dependencies or sync issues here.
  for (m = 0; m < 5; m++) {
    sum[m] = 0.0;
  }

  // 1. Data Dependencies & Synchronization Overhead (cont.):
  // The main computation loop performs accumulation into sum[m].
  // This is a reduction operation. The original code used local temporaries
  // (sum0..sum4) and added them at the end, which is a form of manual privatization
  // and reduction. To better align with OpenMP's efficient 'reduction' clause,
  // we modify the loop body to accumulate directly into the 'sum' array elements.
  // This structure is directly suitable for OpenMP reduction(+:sum).
  // 2. Load Imbalance:
  // The workload within the nested loops is uniform per (i, j, k) iteration.
  // Parallelizing either the 'i' loop (from ist to iend) or the 'j' loop
  // (from jst to jend) allows for distributing the work across threads.
  // The current nested loop structure is well-suited for applying OpenMP parallel for
  // on the 'i' or 'j' loop.
  for (i = ist; i <= iend; i++) {
    for (j = jst; j <= jend; j++) {
      for (k = 1; k <= nz0-2; k++) {
        // Accumulate the squares of elements of v into the sum array.
        // This pattern (sum[m] = sum[m] + term) is recognized and optimized
        // by OpenMP's reduction clause for efficient parallel accumulation.
	sum[0] = sum[0] + v[i][j][k][0] * v[i][j][k][0];
	sum[1] = sum[1] + v[i][j][k][1] * v[i][j][k][1];
	sum[2] = sum[2] + v[i][j][k][2] * v[i][j][k][2];
	sum[3] = sum[3] + v[i][j][k][3] * v[i][j][k][3];
	sum[4] = sum[4] + v[i][j][k][4] * v[i][j][k][4];
      }
    }
  }

  // Calculate the normalization factor.
  // This calculation is independent of the loops above and should be done once.
  norm_factor = (double)(nx0-2)*(ny0-2)*(nz0-2);

  // 1. Data Dependencies & Synchronization Overhead (cont.):
  // The final loop calculates the L2 norm components. Each iteration (for a given m)
  // is independent of other iterations.
  // 2. Load Imbalance: This loop is very small (5 iterations), so imbalance is minimal.
  // 3. Synchronization Overhead: No synchronization is required within this loop
  // as iterations are independent and operate on elements of 'sum' that are
  // assumed to be finalized by the previous computation step. This loop can be
  // parallelized independently after the main reduction loop completes.
  for (m = 0;  m < 5; m++) {
    // Perform the division and square root.
    // Add a check to avoid division by zero or taking the sqrt of a negative number
    // (though sum[m] should be non-negative from the squares).
    if (norm_factor > 0 && sum[m] >= 0.0) {
        sum[m] = sqrt ( sum[m] / norm_factor );
    } else {
        // Handle edge cases, e.g., if the problem size results in a zero norm_factor.
        // Setting to 0.0 is a common approach, depending on the specific requirements.
        sum[m] = 0.0;
    }
  }
}