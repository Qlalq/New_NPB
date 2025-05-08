#include <stdio.h>

static void binvrhs_refined(double lhs[5][5], double r[5]) {
  // This function performs steps similar to Gaussian elimination
  // to transform the matrix lhs and vector r.

  // The outer loop iterates through the pivot steps (k = 0 to 4).
  // This loop is inherently sequential because the transformations
  // for step k depend on the results of step k-1.
  for (int k = 0; k < 5; ++k) {

    // Step 1: Normalize the pivot row 'k'.
    // This makes the diagonal element lhs[k][k] equal to 1.
    // The inverse of the pivot element is calculated once.
    // Assuming lhs[k][k] is non-zero.
    double pivot = 1.00 / lhs[k][k];

    // Normalize elements in row k from column k+1 to 4.
    // These operations within the row are independent of each other.
    // This loop could potentially be parallelized, but is very small for N=5.
    for (int j = k + 1; j < 5; ++j) {
      lhs[k][j] = lhs[k][j] * pivot; // Using assignment operator explicitly as in original
    }

    // Normalize the corresponding element in the right-hand side vector r.
    r[k] = r[k] * pivot; // Using assignment operator explicitly as in original

    // Step 2: Eliminate elements in column k for all other rows 'i'.
    // This step uses the normalized pivot row 'k' to make the elements
    // lhs[i][k] zero for all i != k.
    // The operations for different rows 'i' are independent of each other
    // for a fixed pivot step 'k'. This loop is the main target for parallelization
    // using OpenMP's #pragma omp parallel for.
    for (int i = 0; i < 5; ++i) {
      // Skip the pivot row itself (no elimination needed for lhs[k][k] as it's 1).
      if (i != k) {
        // The coefficient used to eliminate lhs[i][k] in row i.
        // This coefficient is calculated once per row i based on the value
        // of lhs[i][k] before this pivot step k.
        double coeff = lhs[i][k];

        // Update elements in row i from column k+1 to 4.
        // Subtract a multiple of the normalized pivot row k from row i.
        // These operations within row i are independent for different j.
        // Parallelizing the outer loop over i is generally preferred
        // for load balancing and reduced overhead for N=5.
        for (int j = k + 1; j < 5; ++j) {
          lhs[i][j] = lhs[i][j] - coeff * lhs[k][j]; // Using assignment operator explicitly as in original
        }
        // Update the corresponding element in the right-hand side vector r.
        r[i] = r[i] - coeff * r[k]; // Using assignment operator explicitly as in original
      }
    }
  }
  // After the main loop completes for k=0 to 4, the matrix lhs is transformed,
  // and the vector r contains the solution or a transformed vector.
}