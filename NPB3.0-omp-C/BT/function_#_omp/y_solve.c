static void y_solve(void) {
  lhsy();
  y_solve_cell();
  y_backsubstitute();
}