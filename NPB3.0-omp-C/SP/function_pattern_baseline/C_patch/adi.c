static void adi(void) {
 #pragma omp master
 {
   compute_rhs();
   txinvr();
   x_solve();
   y_solve();
   z_solve();
   add();
 }
}