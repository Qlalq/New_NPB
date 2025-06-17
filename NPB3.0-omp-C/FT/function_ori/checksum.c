static void checksum(int i, dcomplex u1[NZ][NY][NX], int d[3]) {

#pragma omp parallel default(shared) 
{

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

    int j, q,r,s, ierr;
    dcomplex chk,allchk;
    
    chk.real = 0.0;
    chk.imag = 0.0;


#pragma omp for nowait
    for (j = 1; j <= 1024; j++) {
	q = j%NX+1;
	if (q >= xstart[0] && q <= xend[0]) {
            r = (3*j)%NY+1;
            if (r >= ystart[0] && r <= yend[0]) {
		s = (5*j)%NZ+1;
		if (s >= zstart[0] && s <= zend[0]) {
		  cadd(chk,chk,u1[s-zstart[0]][r-ystart[0]][q-xstart[0]]);
		}
	    }
	}
    }
#pragma omp critical
    {
	sums[i].real += chk.real;
	sums[i].imag += chk.imag;
    }
#pragma omp barrier
#pragma omp single
  {    
    /* complex % real */
    sums[i].real = sums[i].real/(double)(NTOTAL);
    sums[i].imag = sums[i].imag/(double)(NTOTAL);

    printf("T = %5d     Checksum = %22.12e %22.12e\n",
	   i, sums[i].real, sums[i].imag);
  }
}
}