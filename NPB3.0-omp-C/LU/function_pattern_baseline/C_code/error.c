static void error(void) {
    int i, j, k, m;
    int iglob, jglob;
    double tmp;
    double u000ijk[5];

    /* initialize global error norms */
    for (m = 0; m < 5; m++) {
        errnm[m] = 0.0;
    }

    /* parallelize over i and j with reduction on errnm[0..4] */
#pragma omp parallel for collapse(2) schedule(static) \
    private(iglob, jglob, tmp, u000ijk) \
    reduction(+: errnm[0], errnm[1], errnm[2], errnm[3], errnm[4])
    for (i = ist; i <= iend; i++) {
        for (j = jst; j <= jend; j++) {
            iglob = i;
            jglob = j;
            for (k = 1; k <= nz-2; k++) {
                exact(iglob, jglob, k, u000ijk);
                for (m = 0; m < 5; m++) {
                    tmp = u000ijk[m] - u[i][j][k][m];
                    errnm[m] += tmp * tmp;
                }
            }
        }
    }

    /* finalize root-mean-square norms */
    for (m = 0; m < 5; m++) {
        errnm[m] = sqrt(errnm[m] / ((nx0-2)*(ny0-2)*(nz0-2)));
    }
}