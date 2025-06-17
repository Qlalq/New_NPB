static void fftz2 (int is, int l, int m, int n, int ny, int ny1,
		   dcomplex u[NX], dcomplex x[NX][FFTBLOCKPAD],
		   dcomplex y[NX][FFTBLOCKPAD]) {
    int k,n1,li,lj,lk,ku,i,j,i11,i12,i21,i22;
    dcomplex u1,x11,x21;
    n1 = n / 2;
    if (l-1 == 0) {
	lk = 1;
    } else {
	lk = 2 << ((l - 1)-1);
    }
    if (m-l == 0) {
	li = 1;
    } else {
	li = 2 << ((m - l)-1);
    }
    lj = 2 * lk;
    ku = li;
#pragma omp parallel for private(i11, i12, i21, i22, u1)
    for (i = 0; i < li; i++) {
        i11 = i * lk;
        i12 = i11 + n1;
        i21 = i * lj;
        i22 = i21 + lk;
        if (is >= 1) {
          u1.real = u[ku+i].real;
          u1.imag = u[ku+i].imag;
        } else {
          u1.real = u[ku+i].real;
          u1.imag = -u[ku+i].imag;
        }
        for (k = 0; k < lk; k++) {
	    for (j = 0; j < ny; j++) {
		double x11real, x11imag;
		double x21real, x21imag;
		x11real = x[i11+k][j].real;
		x11imag = x[i11+k][j].imag;
		x21real = x[i12+k][j].real;
		x21imag = x[i12+k][j].imag;
		y[i21+k][j].real = x11real + x21real;
		y[i21+k][j].imag = x11imag + x21imag;
		y[i22+k][j].real = u1.real * (x11real - x21real)
		    - u1.imag * (x11imag - x21imag);
		y[i22+k][j].imag = u1.real * (x11imag - x21imag)
		    + u1.imag * (x11real - x21real);
	    }
	}
    }
}