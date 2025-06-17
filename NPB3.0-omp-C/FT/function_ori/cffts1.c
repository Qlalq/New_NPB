static void cffts1(int is, int d[3], dcomplex x[NZ][NY][NX],
		   dcomplex xout[NZ][NY][NX],
		   dcomplex y0[NX][FFTBLOCKPAD],
		   dcomplex y1[NX][FFTBLOCKPAD]) {

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

    int logd[3];
    int i, j, k, jj;

    for (i = 0; i < 3; i++) {
	logd[i] = ilog2(d[i]);
    }

#pragma omp parallel default(shared) private(i,j,k,jj) shared(is)
{
dcomplex y0[NX][FFTBLOCKPAD];
dcomplex y1[NX][FFTBLOCKPAD];

#pragma omp for 	
    for (k = 0; k < d[2]; k++) {
	for (jj = 0; jj <= d[1] - fftblock; jj+=fftblock) {
/*          if (TIMERS_ENABLED == TRUE) timer_start(T_FFTCOPY); */
            for (j = 0; j < fftblock; j++) {
		for (i = 0; i < d[0]; i++) {

		    y0[i][j].real = x[k][j+jj][i].real;

		    y0[i][j].imag = x[k][j+jj][i].imag;
		}
	    }
/*          if (TIMERS_ENABLED == TRUE) timer_stop(T_FFTCOPY); */
            
/*          if (TIMERS_ENABLED == TRUE) timer_start(T_FFTLOW); */
            cfftz (is, logd[0],
		   d[0], y0, y1);
	    
/*          if (TIMERS_ENABLED == TRUE) timer_stop(T_FFTLOW); */
/*          if (TIMERS_ENABLED == TRUE) timer_start(T_FFTCOPY); */
            for (j = 0; j < fftblock; j++) {
		for (i = 0; i < d[0]; i++) {
		  xout[k][j+jj][i].real = y0[i][j].real;
		  xout[k][j+jj][i].imag = y0[i][j].imag;
		}
	    }
/*          if (TIMERS_ENABLED == TRUE) timer_stop(T_FFTCOPY); */
	}
    }
}
}