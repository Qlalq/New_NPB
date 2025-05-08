static void fft_init (int n) {
    int m,nu,ku,i,j,ln;
    double t, ti;
    nu = n;
    m = ilog2(n);
    u[0].real = (double)m;
    u[0].imag = 0.0;
    ku = 1;
    ln = 1;
    for (j = 1; j <= m; j++) {
	t = PI / ln;
#pragma omp parallel for default(shared) private(i, ti)
	for (i = 0; i <= ln - 1; i++) {
            ti = i * t;
            u[i+ku].real = cos(ti);
	    u[i+ku].imag = sin(ti);
	}
	ku = ku + ln;
	ln = 2 * ln;
    }
}