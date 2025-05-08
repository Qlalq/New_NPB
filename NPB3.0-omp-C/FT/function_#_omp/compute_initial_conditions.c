static void compute_initial_conditions(dcomplex u0[NZ][NY][NX], int d[3]) {
    int k;
    double x0, start, an, dummy;
    static double tmp[NX*2*MAXDIM+1];
    int i, j, t;
    start = SEED;
    ipow46(A, (zstart[0] - 1) * 2 * NX * NY + (ystart[0] - 1) * 2 * NX, &an);
    dummy = randlc(&start, an);
    ipow46(A, 2 * NX * NY, &an);
    {
        for (k = 0; k < dims[0][2]; k++) {
            x0 = start;
            vranlc(2 * NX * dims[0][1], &x0, A, tmp);
            t = 1;
            for (j = 0; j < dims[0][1]; j++) {
                for (i = 0; i < NX; i++) {
                    u0[k][j][i].real = tmp[t++];
                    u0[k][j][i].imag = tmp[t++];
                }
            }
            if (k != dims[0][2]) {
                {
                    dummy = randlc(&start, an);
                }
            }
        }
    }
}