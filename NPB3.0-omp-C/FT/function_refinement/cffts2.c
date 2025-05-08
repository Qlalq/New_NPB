static void cffts2(int is, int d[3], dcomplex x[NZ][NY][NX],
		   dcomplex xout[NZ][NY][NX],
		   dcomplex y0_param[NX][FFTBLOCKPAD], // Parameter kept as per instruction, but will be unused internally
		   dcomplex y1_param[NX][FFTBLOCKPAD]) { // Parameter kept as per instruction, but will be unused internally
    int logd[3];
    int i, j, k, ii;

    // Calculate logd values. This loop has no dependencies and can be run sequentially.
    for (i = 0; i < 3; i++) {
	logd[i] = ilog2(d[i]);
    }

    // The original code structure processes distinct blocks of the 3D data.
    // These blocks are defined by the outer loop index 'k' (across d[2])
    // and the block starting index 'ii' (across d[0] in steps of fftblock).
    // Each (k, ii) block processing is largely independent, especially regarding
    // the temporary buffers 'y0' and 'y1' used within the block via the cfftz call.

    // To prepare for efficient OpenMP parallelization, we need to ensure that
    // each parallel task processing a (k, ii) block has its own private copies
    // of the temporary buffers (y0 and y1) used for the FFT computation.
    // This decouples data dependencies between parallel tasks.

    // The loops over k and ii are the primary candidates for parallelization.
    // We will declare the temporary buffers 'y0_local' and 'y1_local' inside
    // one of these loops (e.g., the 'ii' loop) to ensure they are private
    // for each iteration (and thus for each parallel task assigned an iteration).

    for (k = 0; k < d[2]; k++) {
        // This loop iterates over the third dimension (k).
        // It is a good candidate for parallelization.
        for (ii = 0; ii <= d[0] - fftblock; ii+=fftblock) {
            // This loop iterates over blocks along the first dimension (ii).
            // It is also a good candidate for parallelization.
            // Processing within each (k, ii) iteration is self-contained,
            // operating on specific input/output data blocks and temporary buffers.

            // Declare local, private temporary buffers for this (k, ii) block.
            // These buffers replace the shared parameters y0_param and y1_param
            // to eliminate data hazards when k or ii loops are parallelized.
            // The size matches the dimensions required by cfftz and the inner loops.
            dcomplex y0_local[NX][FFTBLOCKPAD]; // Private temporary buffer for this task
            dcomplex y1_local[NX][FFTBLOCKPAD]; // Private temporary buffer for this task

	    // Copy data from the relevant block in the input array 'x'
            // into the private temporary buffer 'y0_local'.
            // This copy loop operates on the current (k, ii) block data only.
	    for (j = 0; j < d[1]; j++) { // Loop over the second dimension (d[1])
		for (i = 0; i < fftblock; i++) { // Loop over the block size in the first dimension (fftblock)
		    y0_local[j][i].real = x[k][j][i+ii].real;
		    y0_local[j][i].imag = x[k][j][i+ii].imag;
		}
	    }

	    // Perform the core FFT computation on the private temporary buffer.
	    // This call operates on y0_local and uses y1_local as scratch space.
	    // The FFT is along the dimension corresponding to d[1] (index j).
	    cfftz (is, logd[1], // log of the dimension size
		   d[1],        // dimension size
                   y0_local,    // Input/Output buffer (private)
                   y1_local);   // Temporary scratch buffer (private)

            // Copy the results from the private temporary buffer 'y0_local'
            // back into the relevant block in the output array 'xout'.
            // This copy loop operates only on the current (k, ii) block data.
           for (j = 0; j < d[1]; j++) { // Loop over the second dimension (d[1])
	       for (i = 0; i < fftblock; i++) { // Loop over the block size in the first dimension (fftblock)
		   xout[k][j][i+ii].real = y0_local[j][i].real;
		   xout[k][j][i+ii].imag = y0_local[j][i].imag;
	       }
	   }
	} // End of loop over ii blocks
    } // End of loop over k slices

    // The parameters y0_param and y1_param from the function signature
    // are effectively unused in this refactored version. This approach of
    // using local, private buffers (y0_local, y1_local) for each task
    // (each iteration of the parallelized loop) ensures that there are no
    // data dependencies on shared temporary space, thereby addressing:
    // 1. Data Dependencies: Decoupled by privatization of y0/y1.
    // 2. Load Imbalance: Workload per (k, ii) block is relatively uniform,
    //    suitable for static or dynamic scheduling of k or ii loops.
    // 3. Synchronization Overhead: Minimized by avoiding shared mutable state (y0/y1).
}