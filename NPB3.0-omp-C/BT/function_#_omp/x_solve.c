static void x_solve(void) {
  lhsx();
  x_solve_cell();
  x_backsubstitute();
}