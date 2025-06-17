static void compute_indexmap(int indexmap[NZ][NY][NX], int d[3]) {

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c compute function from local (i,j,k) to ibar^2+jbar^2+kbar^2 
c for time evolution exponent. 
c-------------------------------------------------------------------*/

    int i, j, k, ii, ii2, jj, ij2, kk;
    double ap;

/*--------------------------------------------------------------------
c basically we want to convert the fortran indices 
c   1 2 3 4 5 6 7 8 
c to 
c   0 1 2 3 -4 -3 -2 -1
c The following magic formula does the trick:
c mod(i-1+n/2, n) - n/2
c-------------------------------------------------------------------*/

#pragma omp parallel for default(shared) private(i,j,k,ii,ii2,jj,ij2,kk)    
    for (i = 0; i < dims[2][0]; i++) {
	ii =  (i+1+xstart[2]-2+NX/2)%NX - NX/2;
	ii2 = ii*ii;
	for (j = 0; j < dims[2][1]; j++) {
            jj = (j+1+ystart[2]-2+NY/2)%NY - NY/2;
            ij2 = jj*jj+ii2;
            for (k = 0; k < dims[2][2]; k++) {
		kk = (k+1+zstart[2]-2+NZ/2)%NZ - NZ/2;
		indexmap[k][j][i] = kk*kk+ij2;
	    }
	}
    }

/*--------------------------------------------------------------------
c compute array of exponentials for time evolution. 
c-------------------------------------------------------------------*/
    ap = - 4.0 * ALPHA * PI * PI;

    ex[0] = 1.0;
    ex[1] = exp(ap);
    for (i = 2; i <= EXPMAX; i++) {
	ex[i] = ex[i-1]*ex[1];
    }

}