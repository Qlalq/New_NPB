static void adi(void) {
 #pragma omp single
 {
    compute_rhs();
    x_solve();
    y_solve();
    z_solve();
    add();
 }
}