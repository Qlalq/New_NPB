static void cffts3(int is, int d[3], dcomplex x[NZ][NY][NX],
		   dcomplex xout[NZ][NY][NX],
		   dcomplex y0[NX][FFTBLOCKPAD],
		   dcomplex y1[NX][FFTBLOCKPAD]) {
    int logd[3];
    int i, j, k, ii;
    for (i = 0;i < 3; i++) {
	logd[i] = ilog2(d[i]);
    }
{
dcomplex y0[NX][FFTBLOCKPAD];
dcomplex y1[NX][FFTBLOCKPAD];
    for (j = 0; j < d[1]; j++) {
        for (ii = 0; ii <= d[0] - fftblock; ii+=fftblock) {
	    for (k = 0; k < d[2]; k++) {
		for (i = 0; i < fftblock; i++) {
		    y0[k][i].real = x[k][j][i+ii].real;
		    y0[k][i].imag = x[k][j][i+ii].imag;
		}
	    }
           cfftz (is, logd[2],
		  d[2], y0, y1);
           for (k = 0; k < d[2]; k++) {
	       for (i = 0; i < fftblock; i++) {
		   xout[k][j][i+ii].real = y0[k][i].real;
		   xout[k][j][i+ii].imag = y0[k][i].imag;
	       }
	   }
	}
    }
}
}