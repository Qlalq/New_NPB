static void cffts1(int is, int d[3], dcomplex x[NZ][NY][NX],
		   dcomplex xout[NZ][NY][NX])
{
    int logd[3];
    int i, j, k, jj;

    // Calculate logd for each dimension
    for (i = 0; i < 3; i++) {
	logd[i] = ilog2(d[i]);
    }

    // Declare temporary arrays y0 and y1 here.
    // These arrays are temporary buffers used within the nested loops below.
    // For parallel execution of the outer loops (k or jj), these arrays
    // must be private to each thread's work unit to avoid data races.
    // Declaring them here at the function scope, before the loops,
    // makes their scope clear and prepares them for privatization
    // using OpenMP clauses (e.g., #pragma omp for private(y0, y1)).
    // Note: The original code had y0 and y1 as parameters that were shadowed
    // by local declarations inside a block. This version removes the unused
    // parameters and the confusing block, declaring the required temporaries locally.
    dcomplex y0[NX][FFTBLOCKPAD];
    dcomplex y1[NX][FFTBLOCKPAD];


    // The outermost loop iterates through the third dimension (d[2]).
    // Iterations of this loop are independent regarding the input x and output xout,
    // as they operate on distinct planes (k). This loop is a primary candidate
    // for parallelization.
    for (k = 0; k < d[2]; k++) {

        // The inner loop iterates through blocks in the second dimension (d[1]).
        // Iterations of this loop are also independent for a fixed k plane,
        // processing different blocks of data within that plane. This loop
        // can also be parallelized, either instead of or in addition to the k loop,
        // provided y0 and y1 are handled correctly (e.g., privatized per jj iteration).
        for (jj = 0; jj <= d[1] - fftblock; jj+=fftblock) {

            // Copy data from x (input) to y0 (temporary buffer) for the current block.
            // These inner loops perform structured data movement.
            for (j = 0; j < fftblock; j++) {
		        for (i = 0; i < d[0]; i++) {
		            y0[i][j].real = x[k][j+jj][i].real;
		            y0[i][j].imag = x[k][j+jj][i].imag;
		        }
	        }

            // Perform the 1D FFT computation on the temporary buffer y0.
            // The result is typically written back to y0 or into y1 depending on cfftz.
            // This function call represents the core computation on the block data.
            cfftz (is, logd[0],
		           d[0], y0, y1);

            // Copy data from y0 (temporary buffer containing FFT result)
            // to xout (output) for the current block.
            // These inner loops perform structured data movement.
            for (j = 0; j < fftblock; j++) {
		        for (i = 0; i < d[0]; i++) {
		            xout[k][j+jj][i].real = y0[i][j].real; // Assuming cfftz writes result to y0
		            xout[k][j+jj][i].imag = y0[i][j].imag; // Assuming cfftz writes result to y0
		        }
	        }
	    }
    }
}