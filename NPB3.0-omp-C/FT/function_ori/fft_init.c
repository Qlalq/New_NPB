static void fft_init (int n) {

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c compute the roots-of-unity array that will be used for subsequent FFTs. 
c-------------------------------------------------------------------*/

    int m,nu,ku,i,j,ln;
    double t, ti;


/*--------------------------------------------------------------------
c   Initialize the U array with sines and cosines in a manner that permits
c   stride one access at each FFT iteration.
c-------------------------------------------------------------------*/
    nu = n;
    m = ilog2(n);
    u[0].real = (double)m;
    u[0].imag = 0.0;
    ku = 1;
    ln = 1;

    for (j = 1; j <= m; j++) {
	t = PI / ln;
         
	for (i = 0; i <= ln - 1; i++) {
            ti = i * t;
            u[i+ku].real = cos(ti);
	    u[i+ku].imag = sin(ti);
	}
         
	ku = ku + ln;
	ln = 2 * ln;
    }
}