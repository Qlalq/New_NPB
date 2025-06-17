#include <math.h>
#include <omp.h>

static void l2norm (int nx0, int ny0, int nz0,
                    int ist, int iend,
                    int jst, int jend,
                    double v[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5],
                    double sum[5]) {
    int i, j, k, m;
    double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;

    /* initialize output */
    for (m = 0; m < 5; m++) {
        sum[m] = 0.0;
    }

    /* parallel region: accumulate squares into per-thread sum0…sum4, then reduce */
    #pragma omp parallel reduction(+:sum0,sum1,sum2,sum3,sum4)
    {
        #pragma omp for collapse(3) schedule(static)
        for (i = ist; i <= iend; i++) {
            for (j = jst; j <= jend; j++) {
                for (k = 1; k <= nz0-2; k++) {
                    sum0 += v[i][j][k][0] * v[i][j][k][0];
                    sum1 += v[i][j][k][1] * v[i][j][k][1];
                    sum2 += v[i][j][k][2] * v[i][j][k][2];
                    sum3 += v[i][j][k][3] * v[i][j][k][3];
                    sum4 += v[i][j][k][4] * v[i][j][k][4];
                }
            }
        }
    }  /* implicit barrier and reduction here */

    /* accumulate into sum array */
    sum[0] += sum0;
    sum[1] += sum1;
    sum[2] += sum2;
    sum[3] += sum3;
    sum[4] += sum4;

    /* normalize and take square root */
    {
        double denom = 1.0 / ((double)(nx0-2) * (ny0-2) * (nz0-2));
        for (m = 0; m < 5; m++) {
            sum[m] = sqrt(sum[m] * denom);
        }
    }
}