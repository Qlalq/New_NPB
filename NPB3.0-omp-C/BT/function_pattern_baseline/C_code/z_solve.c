static void z_solve(void) {
    #pragma omp parallel sections
    {
        #pragma omp section
        lhsz();

        #pragma omp section
        z_solve_cell();

        #pragma omp section
        z_backsubstitute();
    }
}