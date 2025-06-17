static void x_solve(void) {
  // Create a new parallel region.
  // This allows lhsx(), x_solve_cell(), and x_backsubstitute()
  // to use worksharing constructs like #pragma omp for if they
  // are designed to do so, utilizing the threads of this new team.
  // This assumes nested parallelism is enabled and beneficial.
  #pragma omp parallel
  {
    // Ensure that the sequence of calls is executed by a single thread
    // within the newly created parallel team.
    // #pragma omp single ensures the block is executed by one arbitrary thread
    // from the team and includes an implicit barrier at the end (unless 'nowait' is specified).
    // #pragma omp master could also be used; it's executed by the master thread
    // and has no implicit barrier. Given the functions are black boxes, 'single'
    // is a safer default choice.
    #pragma omp single
    {
      lhsx();
      x_solve_cell();
      x_backsubstitute();
    }
  } // Implicit barrier: all threads in the new team synchronize here.
}