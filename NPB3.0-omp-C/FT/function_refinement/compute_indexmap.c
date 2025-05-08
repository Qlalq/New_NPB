#include <stdio.h>
#include <math.h>

// Assuming these variables/constants are defined elsewhere and accessible
// extern int dims[3][3]; // e.g., dims[2][0], dims[2][1], dims[2][2]
// extern int xstart[3], ystart[3], zstart[3]; // e.g., xstart[2], ystart[2], zstart[2]
// extern int NX, NY, NZ, EXPMAX; // Dimensions and loop bound
// extern double ALPHA, PI; // Constants
// extern double ex[]; // Array to fill, size EXPMAX + 1

static void compute_indexmap_refined(int indexmap[NZ][NY][NX], int d[3]) {
    int i, j, k, ii, ii2, jj, ij2, kk;
    double ap;

    // Part 1: Compute indexmap
    // This nested loop is already structured for parallelism.
    // Each iteration (i, j, k) computes a value for a unique indexmap element indexmap[k][j][i].
    // There are no loop-carried dependencies within this nested loop structure.
    // Load balancing is generally good if the dimensions dims[2][0], dims[2][1], dims[2][2] are reasonable.
    // This part is ready for direct OpenMP parallelization (e.g., using collapse).
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

    // Part 2: Compute exponential array ex
    ap = - 4.0 * ALPHA * PI * PI;

    ex[0] = 1.0;
    // The original loop was:
    // for (i = 2; i <= EXPMAX; i++) { ex[i] = ex[i-1]*ex[1]; }
    // This loop has a loop-carried dependency: ex[i] depends on ex[i-1].
    // The recurrence relation ex[i] = ex[i-1] * ex[1] with ex[1] = exp(ap)
    // results in ex[i] = ex[1]^i = (exp(ap))^i = exp(ap * i) for i >= 1.
    // We can compute each ex[i] independently using the exp function from math.h.
    // This eliminates the loop-carried dependency, making the loop fully parallelizable.
    // This also avoids potential synchronization issues if parallelized naively.
    for (i = 1; i <= EXPMAX; i++) {
        ex[i] = exp(ap * (double)i);
    }
}