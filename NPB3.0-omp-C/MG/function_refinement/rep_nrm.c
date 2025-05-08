#include <stdio.h>

// Assuming nx, ny, nz are globally accessible arrays as in the original context
extern int nx[], ny[], nz[];

// Assuming norm2u3 is an existing function with the following signature:
// void norm2u3(double ***u, int n1, int n2, int n3, double *rnm2, double *rnmu, int nx, int ny, int nz);
// (Declaration might be needed if not in a header file)
extern void norm2u3(double ***u, int n1, int n2, int n3, double *rnm2, double *rnmu, int nx_val, int ny_val, int nz_val);


// Define a struct to explicitly hold the results of the norm computation.
// This helps to structure the data flow, making the dependency of the
// reporting phase on the computation results explicit.
typedef struct {
    double rnm2;
    double rnmu;
} ComputedNorms;


static void rep_nrm(double ***u, int n1, int n2, int n3,
		    char *title, int kk) {

    // Use a local struct variable to store the results of the computation phase.
    // This variable is private to this function call instance.
    ComputedNorms norms;

    // --- Phase 1: Computation of norms ---
    // This call performs the main computational work.
    // The results are written into the 'norms' struct members via pointers.
    // Potential OpenMP parallelization could be applied *within* the
    // 'norm2u3' function implementation itself, or by parallelizing an
    // outer loop that calls 'rep_nrm' for different 'kk' values.
    // No loop-carried dependencies within this specific function's body
    // in this phase.
    norm2u3(u, n1, n2, n3, &norms.rnm2, &norms.rnmu, nx[kk], ny[kk], nz[kk]);

    // --- Phase 2: Reporting / Output ---
    // This phase uses the computed results stored in the 'norms' struct.
    // The printf call writes to a shared resource (stdout).
    // If this function ('rep_nrm') is executed concurrently by multiple
    // OpenMP threads (e.g., in a parallel loop over 'kk'), this printf
    // call will require synchronization (like an OpenMP 'critical' section
    // or ensuring only one thread prints) to prevent mixed output from
    // different threads.
    // The data being printed (norms.rnm2, norms.rnmu) is local to this
    // function instance and does not cause data races itself, but the
    // act of writing to stdout is a synchronization point.
    printf(" Level%2d in %8s: norms =%21.14e%21.14e\n",
	   kk, title, norms.rnm2, norms.rnmu);

    // This structured approach separates the computation phase (producing
    // results stored in 'norms') from the reporting phase (consuming
    // results from 'norms' and performing I/O), making it clearer where
    // computation parallelism might occur (Phase 1) and where synchronization
    // is needed due to shared resources (Phase 2).
}