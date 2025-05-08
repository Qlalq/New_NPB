#include <stdio.h> // Assuming this is needed based on typical FFT usage
#include <complex.h> // Assuming dcomplex is complex double

// Define dcomplex if not provided by a header (common in older codes or specific frameworks)
#ifndef dcomplex
typedef double complex dcomplex;
#endif

// Assume NX and FFTBLOCKPAD are defined elsewhere (e.g., in a header file)
// #define NX ...
// #define FFTBLOCKPAD ...

static void fftz2_refined (int is, int l, int m, int n, int ny, int ny1,
		   dcomplex u[NX], dcomplex x[NX][FFTBLOCKPAD],
		   dcomplex y[NX][FFTBLOCKPAD]) {
    int k, n1, li, lj, lk, ku, i, j;

    n1 = n / 2;

    // Calculate lk = 2^(l-1)
    if (l-1 == 0) {
	    lk = 1;
    } else {
	    lk = 2 << ((l - 1)-1);
    }

    // Calculate li = 2^(m-l)
    if (m-l == 0) {
	    li = 1;
    } else {
	    li = 2 << ((m - l)-1);
    }

    lj = 2 * lk;
    ku = li;

    // The outermost loop iterates through blocks that are independent of each other.
    // Each iteration of 'i' processes a distinct set of input data
    // (x[i*lk..i*lk+lk-1][:] and x[i*lk+n/2..i*lk+n/2+lk-1][:])
    // and writes to a distinct set of output data
    // (y[i*2*lk..i*2*lk+lk-1][:] and y[(2*i+1)*lk..(2*i+1)*lk+lk-1][:])
    // This structure allows for easy parallelization of the 'i' loop.
    for (i = 0; i < li; i++) {

        // These indices and 'u1' are calculated once per outer loop iteration 'i'
        // and are used only within this iteration's work.
        // They are effectively private to each 'i' iteration.
        int i11 = i * lk;
        int i12 = i11 + n1;
        int i21 = i * lj;
        int i22 = i21 + lk;

        dcomplex u1;
        // Conjugate or not based on 'is' flag
        if (is >= 1) {
            u1 = u[ku+i];
        } else {
            u1 = conj(u[ku+i]); // Using complex.h conj
        }

        // The middle loop iterates through sub-blocks (butterflies) within the current 'i' block.
        // For a fixed 'i', iterations of 'k' process disjoint data elements in 'x' and write
        // to disjoint data elements in 'y'. This loop can also be parallelized,
        // potentially combined with the 'i' loop for larger parallel tasks.
        for (k = 0; k < lk; k++) {

            // The innermost loop iterates through the 'j' dimension (likely orthogonal to the FFT dimension).
            // For fixed 'i' and 'k', iterations of 'j' are independent.
            // They access distinct column elements in 'x' and write to distinct column elements in 'y'.
            for (j = 0; j < ny; j++) {

                // Load input values for the current butterfly
                dcomplex x11 = x[i11+k][j];
                dcomplex x21 = x[i12+k][j];

                // Perform the butterfly calculation
                dcomplex term = u1 * x21; // complex multiplication

                // Store results in the output array
                y[i21+k][j] = x11 + term;
                y[i22+k][j] = x11 - term;
            }
        }
    }
}