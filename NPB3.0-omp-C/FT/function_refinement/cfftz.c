#include <stdio.h>
#include <stdlib.h> // For exit

// Assuming dcomplex type is defined elsewhere
// typedef struct { double real; double imag; } dcomplex;

// Assuming global or file-scope constants for array dimensions
// #define NX ...
// #define FFTBLOCKPAD ...

// Assuming global or file-scope variables used by fftz2 and cfftz
// extern int fftblock;
// extern int fftblockpad;
// extern dcomplex u[]; // Used for initialization check

// Assuming fftz2 function is declared elsewhere
// extern void fftz2 (int is, int l, int m, int n, int fftblock, int fftblockpad, dcomplex u[], dcomplex x[][FFTBLOCKPAD], dcomplex y[][FFTBLOCKPAD]);


// Placeholder definitions for compilation if not provided by environment headers
// In a real scenario, these would come from included header files
#ifndef NX
#define NX 1024 // Example placeholder value
#endif
#ifndef FFTBLOCKPAD
#define FFTBLOCKPAD 1024 // Example placeholder value
#endif

// Define the complex number structure if not provided
typedef struct {
    double real;
    double imag;
} dcomplex;

// Declare external variables used (assuming they are global or file-scope static)
extern int fftblock;
extern int fftblockpad;
extern dcomplex u[]; // Assuming u is an array of dcomplex, size at least 1

// Declare external function fftz2
extern void fftz2 (int is, int l, int m, int n, int fftblock, int fftblockpad, dcomplex u[], dcomplex x[][FFTBLOCKPAD], dcomplex y[][FFTBLOCKPAD]);


static void cfftz (int is, int m, int n, dcomplex x[NX][FFTBLOCKPAD],
		   dcomplex y[NX][FFTBLOCKPAD]) {
    int i,j,l,mx;

    // --- Sequential Part: Parameter validation ---
    // This part performs checks and does not involve significant computation on the main data
    // and has no data dependencies preventing parallelization of *itself*.
    // It must complete before the main computation begins.
    mx = (int)(u[0].real); // Access to global u
    if ((is != 1 && is != -1) || m < 1 || m > mx) {
	printf("CFFTZ: Either U has not been initialized, or else\n"
	       "one of the input parameters is invalid%5d%5d%5d\n",
	       is, m, mx);
	exit(1);
    }

    // --- Sequential Part: Main FFT Computation Stages ---
    // This loop processes pairs of FFT stages (l and l+1).
    // The output of fftz2 for stage l becomes the input for fftz2 for stage l+1
    // within the same iteration of this loop.
    // The output of stage l+1 (in one iteration) becomes the input for stage l+2
    // (in the next iteration of the l+=2 loop).
    // This creates a strong loop-carried dependency across iterations (l, l+2, l+4, ...),
    // making the main `for` loop inherently sequential at this level.
    // Data Dependencies: Output buffer from previous stage pair (buf_A) is input for current stage pair.
    // Load Imbalance: Work per iteration (two fftz2 calls) is roughly constant.
    // Synchronization: Dependency prevents simple parallelization; synchronization would be complex and likely negate gains.

    // Use pointers to explicitly manage the two primary buffers (x and y).
    // This clarifies which buffer is being used as input and output for each stage pair.
    dcomplex (*buf_A)[FFTBLOCKPAD] = x; // Buffer A (initially holds input, receives output after even stages)
    dcomplex (*buf_B)[FFTBLOCKPAD] = y; // Buffer B (initially temporary, receives output after odd stages)

    for (l = 1; l <= m; l+=2) {
        // Process layer l (odd stage): Reads from buf_A, writes result to buf_B
        // buf_A holds the input from the previous combined stage (or initial input x).
        // buf_B receives the output of stage l.
        fftz2 (is, l, m, n, fftblock, fftblockpad, u, buf_A, buf_B);

        // Check if this was the last stage (only happens if m is odd).
        // If l == m, we processed the last odd stage, result is in buf_B. No l+1 stage.
        if (l == m) break;

        // Process layer l+1 (even stage): Reads from buf_B, writes result to buf_A
        // buf_B contains the result of layer l, which is the required input for layer l+1.
        // buf_A receives the output of stage l+1.
        fftz2 (is, l + 1, m, n, fftblock, fftblockpad, u, buf_B, buf_A);

        // After this iteration, the result of layer l+1 is in buf_A.
        // For the next pair (l+2, l+3), buf_A is correctly positioned to be the input buffer.
        // The roles of buf_A and buf_B are consistent across iterations of this l+=2 loop:
        // buf_A holds the state after even stages, buf_B holds the state after odd stages.
    }

    // --- Parallelizable Part: Final Data Copy (if needed) ---
    // If m is odd, the last operation was fftz2(..., buf_A, buf_B) for l=m.
    // The final result is currently in buf_B (which points to y).
    // The function is expected to return the result in x (pointed to by buf_A).
    // If m is even, the last operations were for l=m-1 and l+1=m, ending with fftz2(..., buf_B, buf_A).
    // The final result is already in buf_A (which points to x).
    if (m % 2 == 1) {
        // Copy contents from buf_B (y) to buf_A (x)
        // This nested loop involves independent writes to buf_A elements.
        // It is suitable for OpenMP parallelization (e.g., #pragma omp for on the outer j loop).
        // Data Dependencies: None between iterations of j or i loops.
        // Load Imbalance: Fairly balanced copy operation.
        // Synchronization: Not needed between iterations of these loops.
        for (j = 0; j < n; j++) {
            for (i = 0; i < fftblock; i++) {
                buf_A[j][i].real = buf_B[j][i].real;
                buf_A[j][i].imag = buf_B[j][i].imag;
            }
        }
    }
    // If m is even, the result is already in buf_A (x), so no copy is needed.
}