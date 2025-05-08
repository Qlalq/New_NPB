static void setup(void) {
    int ierr, i, j, fstatus;
    printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version"
	   " - FT Benchmark\n\n");
    niter = NITER_DEFAULT;
    printf(" Size                : %3dx%3dx%3d\n", NX, NY, NZ);
    printf(" Iterations          :     %7d\n", niter);
    for (i = 0;i < 3 ; i++) {
	dims[i][0] = NX;
	dims[i][1] = NY;
	dims[i][2] = NZ;
    }
    #pragma omp parallel for // Loop iterates only 3 times, parallelization overhead likely outweighs benefit.
    for (i = 0; i < 3; i++) {
	xstart[i] = 1;
	xend[i]   = NX;
	ystart[i] = 1;
        yend[i]   = NY;
        zstart[i] = 1;
        zend[i]   = NZ;
    }
    fftblock = FFTBLOCK_DEFAULT;
    fftblockpad = FFTBLOCKPAD_DEFAULT;
    if (fftblock != FFTBLOCK_DEFAULT) fftblockpad = fftblock+3;
}