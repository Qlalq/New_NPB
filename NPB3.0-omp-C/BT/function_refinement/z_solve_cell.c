static void z_solve_cell(void) {
  int i, j, k; // i and j will be derived from ij
  int ksize;
  ksize = grid_points[2]-1;

  // Determine the number of points in the i and j dimensions that are processed.
  // Original loops are from 1 to < grid_points[0]-1 and 1 to < grid_points[1]-1.
  // This means i takes values 1, 2, ..., grid_points[0]-2. Number of values = grid_points[0]-2.
  // This means j takes values 1, 2, ..., grid_points[1]-2. Number of values = grid_points[1]-2.
  const int ni = grid_points[0]-2;
  const int nj = grid_points[1]-2;
  const int total_ij = ni * nj;

  // The original computation is structured into three phases based on the k index:
  // k = 0, k = 1..ksize-1, and k = ksize.
  //
  // Dependency Analysis:
  // - The loops over i and j for a *fixed* k are independent. Iterations for different (i, j)
  //   pairs modify different elements of rhs and lhs. These loops are parallelizable.
  // - The loop over k (from 1 to ksize-1) has loop-carried dependencies:
  //   Calculations at k use results from k-1 (rhs[i][j][k-1] and lhs[i][j][k-1][CC]).
  //   This k loop MUST be executed sequentially.
  //
  // Refinement Strategy:
  // - Keep the sequential structure for the k dimension.
  // - Expose the independent work units in the i and j dimensions by flattening the i,j loops
  //   into a single loop over a combined index 'ij'. This makes applying OpenMP parallelization
  //   to the independent work units straightforward.
  // - The mapping from the flattened index 'ij' back to (i, j) is:
  //   i = 1 + (ij / nj)
  //   j = 1 + (ij % nj)
  //   where 1 is the starting index for i and j.

  // Block 1: k = 0 processing
  // All ij iterations are independent and can be parallelized.
  for (int ij = 0; ij < total_ij; ++ij) {
    // Map flattened index ij back to (i, j)
    i = 1 + (ij / nj);
    j = 1 + (ij % nj);

    // Original code for k = 0
    binvcrhs( lhs[i][j][0][BB],
              lhs[i][j][0][CC],
              rhs[i][j][0] );
  }

  // Block 2: k from 1 to ksize-1 processing (sequential k sweep)
  // The loop over k is sequential due to dependencies.
  // Inside the k loop, ij iterations are independent and can be parallelized.
  for (k = 1; k < ksize; k++) {
      // This inner loop over ij is parallelizable for a fixed k.
      for (int ij = 0; ij < total_ij; ++ij) {
          // Map flattened index ij back to (i, j)
          i = 1 + (ij / nj);
          j = 1 + (ij % nj);

          // Original code for 1 <= k < ksize
          matvec_sub(lhs[i][j][k][AA],
                     rhs[i][j][k-1], rhs[i][j][k]);
          matmul_sub(lhs[i][j][k][AA],
                     lhs[i][j][k-1][CC],
                     lhs[i][j][k][BB]);
          binvcrhs( lhs[i][j][k][BB],
                    lhs[i][j][k][CC],
                    rhs[i][j][k] );
      }
  }

  // Block 3: k = ksize processing
  // All ij iterations are independent and can be parallelized.
  for (int ij = 0; ij < total_ij; ++ij) {
    // Map flattened index ij back to (i, j)
    i = 1 + (ij / nj);
    j = 1 + (ij % nj);

    // Original code for k = ksize
    binvrhs( lhs[i][j][ksize][BB],
             rhs[i][j][ksize] );
  }
}