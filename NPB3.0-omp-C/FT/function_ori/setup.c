static void setup(void) {

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

    int ierr, i, j, fstatus;
      
    printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version"
	   " - FT Benchmark\n\n");

    niter = NITER_DEFAULT;

    printf(" Size                : %3dx%3dx%3d\n", NX, NY, NZ);
    printf(" Iterations          :     %7d\n", niter);

/* 1004 format(' Number of processes :     ', i7)
 1005 format(' Processor array     :     ', i3, 'x', i3)
 1006 format(' WARNING: compiled for ', i5, ' processes. ',
     >       ' Will not verify. ')*/

    for (i = 0;i < 3 ; i++) {
	dims[i][0] = NX;
	dims[i][1] = NY;
	dims[i][2] = NZ;
    }


    for (i = 0; i < 3; i++) {
	xstart[i] = 1;
	xend[i]   = NX;
	ystart[i] = 1;
        yend[i]   = NY;
        zstart[i] = 1;
        zend[i]   = NZ;
    }

/*--------------------------------------------------------------------
c Set up info for blocking of ffts and transposes.  This improves
c performance on cache-based systems. Blocking involves
c working on a chunk of the problem at a time, taking chunks
c along the first, second, or third dimension. 
c
c - In cffts1 blocking is on 2nd dimension (with fft on 1st dim)
c - In cffts2/3 blocking is on 1st dimension (with fft on 2nd and 3rd dims)

c Since 1st dim is always in processor, we'll assume it's long enough 
c (default blocking factor is 16 so min size for 1st dim is 16)
c The only case we have to worry about is cffts1 in a 2d decomposition. 
c so the blocking factor should not be larger than the 2nd dimension. 
c-------------------------------------------------------------------*/

    fftblock = FFTBLOCK_DEFAULT;
    fftblockpad = FFTBLOCKPAD_DEFAULT;

    if (fftblock != FFTBLOCK_DEFAULT) fftblockpad = fftblock+3;
}