static void binvcrhs(double lhs[5][5], double c[5][5], double r[5]) {
  int size = 5; // Define size for clarity in loops

  // The main loop iterates through each column (pivot column).
  // This loop is sequential because the pivot element lhs[k][k] and
  // the normalized pivot row (row k) depend on the computations
  // from previous iterations (k-1).
  for (int k = 0; k < size; k++) {
    double pivot, coeff;

    // Step 1: Normalize the pivot row (row k).
    // Divide row k by the pivot element lhs[k][k] to make lhs[k][k] = 1.0.
    // This step must complete before the elimination step using this row.
    pivot = 1.00 / lhs[k][k];

    // Apply normalization to elements of lhs in the pivot row *after* the pivot column.
    // The original code structure implies only these elements are modified in lhs[k].
    for (int j = k + 1; j < size; j++) {
      lhs[k][j] *= pivot;
    }

    // Apply normalization to elements of c in the pivot row (all columns).
    for (int j = 0; j < size; j++) {
      c[k][j] *= pivot;
    }

    // Apply normalization to the element of r in the pivot row.
    r[k] *= pivot;

    // Step 2: Eliminate non-pivot elements in the current column k.
    // For each row i (where i != k), subtract a multiple of row k
    // to make lhs[i][k] = 0.
    // The operations for different rows i (where i != k) are independent
    // and can be executed concurrently for a fixed k. This loop over 'i'
    // is the primary target for OpenMP parallelization.
    for (int i = 0; i < size; i++) {
      // Skip the pivot row itself, as its element in column k is already 1 (implicitly).
      if (i != k) {
        // The coefficient is the element in row i, column k before the current step's elimination.
        // This value needs to be computed based on the current state of lhs[i][k].
        coeff = lhs[i][k];

        // Subtract coeff * row k from row i for lhs (elements after column k).
        // Note that lhs[i][k] becomes 0 after this logical step, but isn't explicitly set to 0 here,
        // matching the structure implied by the original unrolled code.
        for (int j = k + 1; j < size; j++) {
          lhs[i][j] -= coeff * lhs[k][j];
        }

        // Subtract coeff * row k from row i for c (all elements).
        for (int j = 0; j < size; j++) {
          c[i][j] -= coeff * c[k][j];
        }

        // Subtract coeff * row k from row i for r.
        r[i] -= coeff * r[k];
      }
    }
  }
  // After the loops, lhs should ideally be the identity matrix (with diagonal elements explicitly 1 and off-diagonals implicitly 0
  // based on how elements are used/modified in this specific implementation).
  // c should contain the result of A_inverse * C.
  // r should contain the result of A_inverse * R.
}