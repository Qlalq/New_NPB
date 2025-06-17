static void y_solve(void) {
  lhsy();
#pragma omp parallel
  y_solve_cell();
  y_backsubstitute();
}