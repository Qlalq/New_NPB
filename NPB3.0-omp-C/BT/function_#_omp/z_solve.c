static void z_solve(void) {
  lhsz();
  z_solve_cell();
  z_backsubstitute();
}