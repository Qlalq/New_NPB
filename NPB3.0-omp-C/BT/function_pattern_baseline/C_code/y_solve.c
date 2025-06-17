static void y_solve(void) {
  // Create a new parallel region. The single thread that entered y_solve()
  // will become the master of this new team of threads.
  #pragma omp parallel
  {
    // All threads in the team call lhsy().
    // It is assumed that lhsy() contains worksharing constructs (e.g., #pragma omp for)
    // to distribute its workload among the threads of this team.
    // An implicit barrier at the end of the worksharing construct(s) in lhsy()
    // ensures all its work is completed before any thread proceeds.
    lhsy();

    // All threads in the team call y_solve_cell().
    // Similar assumptions about internal worksharing and synchronization apply.
    y_solve_cell();

    // All threads in the team call y_backsubstitute().
    // Similar assumptions about internal worksharing and synchronization apply.
    y_backsubstitute();
  }
  // Implicit barrier at the end of the #pragma omp parallel region.
  // All threads synchronize here, the team is disbanded, and the master thread
  // (the one that originally entered y_solve) continues.
}