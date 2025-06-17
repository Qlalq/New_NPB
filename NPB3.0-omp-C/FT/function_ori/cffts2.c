static void cffts2(int is, int d[3], dcomplex x[NZ][NY][NX],
		   dcomplex xout[NZ][NY][NX],
		   dcomplex y0[NX][FFTBLOCKPAD],
		   dcomplex y1[NX][FFTBLOCKPAD]) {

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

    int logd[3];
    int i, j, k, ii;

    for (i = 0; i < 3; i++) {
	logd[i] = ilog2(d[i]);
    }
#pragma omp parallel default(shared) private(i,j,k,ii) shared(is)
{
dcomplex y0[NX][FFTBLOCKPAD];
dcomplex y1[NX][FFTBLOCKPAD];

#pragma omp for 	
    for (k = 0; k < d[2]; k++) {
        for (ii = 0; ii <= d[0] - fftblock; ii+=fftblock) {
/*	    if (TIMERS_ENABLED == TRUE) timer_start(T_FFTCOPY); */
	    for (j = 0; j < d[1]; j++) {
		for (i = 0; i < fftblock; i++) {
		    y0[j][i].real = x[k][j][i+ii].real;
		    y0[j][i].imag = x[k][j][i+ii].imag;
		}
	    }
/*	    if (TIMERS_ENABLED == TRUE) timer_stop(T_FFTCOPY); */
/*	    if (TIMERS_ENABLED == TRUE) timer_start(T_FFTLOW); */
	    cfftz (is, logd[1], 
		   d[1], y0, y1);
           
/*          if (TIMERS_ENABLED == TRUE) timer_stop(T_FFTLOW); */
/*          if (TIMERS_ENABLED == TRUE) timer_start(T_FFTCOPY); */
           for (j = 0; j < d[1]; j++) {
	       for (i = 0; i < fftblock; i++) {
		   xout[k][j][i+ii].real = y0[j][i].real;
		   xout[k][j][i+ii].imag = y0[j][i].imag;
	       }
	   }
/*           if (TIMERS_ENABLED == TRUE) timer_stop(T_FFTCOPY); */
	}
    }
}
}