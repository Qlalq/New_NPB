static void cfftz (int is, int m, int n, dcomplex x[NX][FFTBLOCKPAD],
		   dcomplex y[NX][FFTBLOCKPAD]) {

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c   Computes NY N-point complex-to-complex FFTs of X using an algorithm due
c   to Swarztrauber.  X is both the input and the output array, while Y is a 
c   scratch array.  It is assumed that N = 2^M.  Before calling CFFTZ to 
c   perform FFTs, the array U must be initialized by calling CFFTZ with IS 
c   set to 0 and M set to MX, where MX is the maximum value of M for any 
c   subsequent call.
c-------------------------------------------------------------------*/

    int i,j,l,mx;

/*--------------------------------------------------------------------
c   Check if input parameters are invalid.
c-------------------------------------------------------------------*/
    mx = (int)(u[0].real);
    if ((is != 1 && is != -1) || m < 1 || m > mx) {
	printf("CFFTZ: Either U has not been initialized, or else\n"
	       "one of the input parameters is invalid%5d%5d%5d\n",
	       is, m, mx);
	exit(1);
    }

/*--------------------------------------------------------------------
c   Perform one variant of the Stockham FFT.
c-------------------------------------------------------------------*/
    for (l = 1; l <= m; l+=2) {
        fftz2 (is, l, m, n, fftblock, fftblockpad, u, x, y);
        if (l == m) break;
	fftz2 (is, l + 1, m, n, fftblock, fftblockpad, u, y, x);
    }

/*--------------------------------------------------------------------
c   Copy Y to X.
c-------------------------------------------------------------------*/
    if (m % 2 == 1) {
	for (j = 0; j < n; j++) {
	    for (i = 0; i < fftblock; i++) {
		x[j][i].real = y[j][i].real;
		x[j][i].imag = y[j][i].imag;
	    }
	}
    }
}