static void z_solve(void) {
  // This function performs a sequence of distinct stages.
  // These stages are likely dependent and must be executed in order.
  // Efficient OpenMP parallelization should focus on identifying and exploiting
  // parallelism *within* each of the following functions where possible,
  // as the calls themselves are sequential due to dependencies.

  // --- Stage 1: Setup Left-Hand Side (LHS) ---
  // This stage typically involves operations to prepare the linear system's
  // left-hand side. Potential for internal data-parallel operations.
  // Dependencies: May depend on global state or input parameters.
  // OpenMP Opportunity: Parallel loops or tasks within lhsz() if independent work units exist.
  lhsz();

  // --- Stage 2: Solve for Cells ---
  // This stage computes intermediate results or solves sub-problems
  // for individual cells or domain partitions.
  // Dependencies: Depends on the setup performed in Stage 1.
  // OpenMP Opportunity: Parallelize over independent cells or partitions within z_solve_cell(). Load balancing may be needed here if cell work varies.
  z_solve_cell();

  // --- Stage 3: Back Substitution / Finalization ---
  // This stage uses the results from Stage 2 to compute final values, often
  // involving backward dependencies (like back-substitution in a matrix solve).
  // Dependencies: Depends on the results from Stage 2. May involve loop-carried dependencies.
  // OpenMP Opportunity: Parallelism within z_backsubstitute() might require careful handling of dependencies (e.g., parallel scan, wavefront) or partitioning if applicable. Minimizing synchronization for write operations is key here.
  z_backsubstitute();
}