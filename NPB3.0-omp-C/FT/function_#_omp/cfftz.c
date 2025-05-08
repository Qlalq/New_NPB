static void cfftz (int is, int m, int n, dcomplex x[NX][FFTBLOCKPAD],
		   dcomplex y[NX][FFTBLOCKPAD]) {
    int i,j,l,mx;
    mx = (int)(u[0].real);
    if ((is != 1 && is != -1) || m < 1 || m > mx) {
	printf("CFFTZ: Either U has not been initialized, or else\n"
	       "one of the input parameters is invalid%5d%5d%5d\n",
	       is, m, mx);
	exit(1);
    }
    for (l = 1; l <= m; l+=2) {
        fftz2 (is, l, m, n, fftblock, fftblockpad, u, x, y);
        if (l == m) break;
	fftz2 (is, l + 1, m, n, fftblock, fftblockpad, u, y, x);
    }
    if (m % 2 == 1) {
	for (j = 0; j < n; j++) {
	    for (i = 0; i < fftblock; i++) {
		x[j][i].real = y[j][i].real;
		x[j][i].imag = y[j][i].imag;
	    }
	}
    }
}