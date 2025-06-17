#include <omp.h>

static void adi(void) {
    #pragma omp parallel sections
    {
        #pragma omp section
        compute_rhs();

        #pragma omp section
        x_solve();

        #pragma omp section
        y_solve();

        #pragma omp section
        z_solve();

        #pragma omp section
        add();
    }
}