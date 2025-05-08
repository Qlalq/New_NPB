static void matvec_sub_refined(double ablock[5][5], double avec[5], double bvec[5]) {
  // Use a temporary array to store the values to be subtracted from bvec.
  // This helps decouple the computation of these values from the modification
  // of bvec, making the independent nature of the computations clearer.
  double subtrahend_vec[5];

  int i;
  // Loop 1: Compute the values to be subtracted for each element of bvec.
  // Each iteration 'i' calculates the dot product of the i-th row of ablock
  // with the vector avec.
  // The computation of subtrahend_vec[i] is completely independent of
  // subtrahend_vec[j] for j != i.
  // This loop is a primary candidate for parallelization across the 'i' dimension.
  for (i = 0; i < 5; i++) {
    double current_row_sub = 0.0; // Accumulator for the dot product for row 'i'
    int j;
    // The inner loop performs a reduction (summation) for a single element
    // of the subtrahend_vec. For small fixed sizes like 5, this is efficient
    // sequentially or can be vectorized by the compiler. For larger sizes,
    // one might consider parallelizing or optimizing this inner loop as well,
    // but the primary gain comes from parallelizing the outer 'i' loop.
    for (j = 0; j < 5; j++) {
      current_row_sub += ablock[i][j] * avec[j];
    }
    // Store the computed value in the temporary array.
    // Each iteration 'i' writes to a unique element subtrahend_vec[i].
    // This write is independent across 'i' iterations.
    subtrahend_vec[i] = current_row_sub;
  }

  // Loop 2: Subtract the pre-computed values from bvec.
  // Each iteration 'i' updates bvec[i] using the value stored in subtrahend_vec[i].
  // The update operation bvec[i] = bvec[i] - subtrahend_vec[i] for a given 'i'
  // is completely independent of the update operation for bvec[j] where j != i.
  // This loop is also highly suitable for parallelization across the 'i' dimension.
  for (i = 0; i < 5; i++) {
    // Read from subtrahend_vec[i] and bvec[i], write to bvec[i].
    // The write operation to bvec[i] is to a unique memory location for each
    // iteration 'i' of this loop. No synchronization is needed for these writes.
    bvec[i] = bvec[i] - subtrahend_vec[i];
  }
}