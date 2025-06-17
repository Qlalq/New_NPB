#include <omp.h>
#include "globals.h"    /* defines nx0, ny0, nz, ce[][13], and u[nz][ny0][nx0][5] */

#pragma omp declare simd
static void exact(int i, int j, int k, double u000ijk[5]) {
    int m;
    double xi, eta, zeta;

    xi   = ((double)i) / (nx0 - 1);
    eta  = ((double)j) / (ny0 - 1);
    zeta = ((double)k) / (nz  - 1);

    /* Vectorize the small m‐loop across the 5 components */
    #pragma omp simd
    for (m = 0; m < 5; m++) {
        u000ijk[m] =
              ce[m][0]
            + ce[m][1] * xi
            + ce[m][2] * eta
            + ce[m][3] * zeta
            + ce[m][4] * xi * xi
            + ce[m][5] * eta * eta
            + ce[m][6] * zeta * zeta
            + ce[m][7] * xi * xi * xi
            + ce[m][8] * eta * eta * eta
            + ce[m][9] * zeta * zeta * zeta
            + ce[m][10] * xi * xi * xi * xi
            + ce[m][11] * eta * eta * eta * eta
            + ce[m][12] * zeta * zeta * zeta * zeta;
    }
}

void compute_exact_solution() {
    int i, j, k;

    /* Replace the three separate parallel‐for pragmas with one 3‐level collapsed loop:
       - collapse(3) lets OpenMP distribute the (k,j,i) iterations in one go
       - schedule(static) gives balanced chunks without runtime cost
       - default(none) forces explicit sharing clauses for safety
    */
    #pragma omp parallel for collapse(3) schedule(static) default(none)              \
        private(i, j, k) shared(u)
    for (k = 0; k < nz; ++k) {
        for (j = 0; j < ny0; ++j) {
            for (i = 0; i < nx0; ++i) {
                exact(i, j, k, u[k][j][i]);
            }
        }
    }
}