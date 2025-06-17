static void adi(void) {
#pragma omp parallel
    compute_rhs();

#pragma omp parallel
    x_solve();

#pragma omp parallel
    y_solve();

#pragma omp parallel
    z_solve();

#pragma omp parallel
    add();
}