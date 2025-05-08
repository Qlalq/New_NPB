static void setup(void) {
    int ierr, i, j, fstatus;

    // Serial initialization of parameters and print statements
    printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version"
	   " - FT Benchmark\n\n");
    niter = NITER_DEFAULT;
    printf(" Size                : %3dx%3dx%3d\n", NX, NY, NZ);
    printf(" Iterations          :     %7d\n", niter);

    // Initialization loop for dims array
    // This loop is suitable for parallelization as each iteration
    // writes to distinct elements and has no loop-carried dependencies.
    for (i = 0;i < 3 ; i++) {
	dims[i][0] = NX;
	dims[i][1] = NY;
	dims[i][2] = NZ;
    }

    // Initialization loop for boundary arrays (xstart, xend, ystart, yend, zstart, zend)
    // This loop is also suitable for parallelization as each iteration
    // writes to distinct elements across multiple arrays and has no
    // loop-carried dependencies.
    for (i = 0; i < 3; i++) {
	xstart[i] = 1;
	xend[i]   = NX;
	ystart[i] = 1;
        yend[i]   = NY;
        zstart[i] = 1;
        zend[i]   = NZ;
    }

    // Serial initialization of FFT parameters
    fftblock = FFTBLOCK_DEFAULT;
    fftblockpad = FFTBLOCKPAD_DEFAULT;
    if (fftblock != FFTBLOCK_DEFAULT) fftblockpad = fftblock+3;
}