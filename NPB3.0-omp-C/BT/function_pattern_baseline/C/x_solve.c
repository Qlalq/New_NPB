static void x_solve(void) {
#pragma omp parallel
#pragma omp single
  {
  lhsx();
  x_solve_cell();
  x_backsubstitute();
  } /* end of single block */
}