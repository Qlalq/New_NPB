#include <stdio.h> // Keep original includes if any

// Assuming the original code intended three separate boundary communication loops
// This refinement separates them and prepares for OpenMP parallelization.
static void comm3_refined(double ***u, int n1, int n2, int n3, int kk) {
    // Although kk is not used in the provided snippet, keep it as per the original signature.

    // This function copies boundary data from inner planes to outer ghost planes.
    // It involves three distinct sets of copies, one for each dimension's faces.
    // The original code had a structural error where the second loop was nested
    // inside the first loop's outermost 'i3' loop. This has been corrected.

    // By placing the three independent loop nests within a single parallel region
    // and using 'nowait', threads can potentially work on different boundary faces
    // concurrently without unnecessary synchronization between the sections.

#pragma omp parallel
    {
        int i1, i2, i3; // Declare loop variables private to each thread

        // Copy faces in the n1 dimension (x-faces: i1=0 and i1=n1-1)
        // These copies are independent for each (i3, i2) pair.
#pragma omp for nowait // Parallelize over i3, no barrier after this loop nest
        for (i3 = 1; i3 < n3 - 1; i3++) {
            for (i2 = 1; i2 < n2 - 1; i2++) {
                u[i3][i2][n1 - 1] = u[i3][i2][1];    // Copy from i1=1 to i1=n1-1
                u[i3][i2][0] = u[i3][i2][n1 - 2];      // Copy from i1=n1-2 to i1=0
            }
        }

        // Copy faces in the n2 dimension (y-faces: i2=0 and i2=n2-1)
        // These copies are independent for each (i3, i1) pair.
        // Assuming n1, n2, n3 >= 3, these write locations (i2=0/n2-1)
        // do not overlap with reads/writes of the first loop (i1=0/n1-1, i1=1/n1-2).
#pragma omp for nowait // Parallelize over i3, no barrier after this loop nest
        for (i3 = 1; i3 < n3 - 1; i3++) {
            for (i1 = 0; i1 < n1; i1++) {
                u[i3][n2 - 1][i1] = u[i3][1][i1];    // Copy from i2=1 to i2=n2-1
                u[i3][0][i1] = u[i3][n2 - 2][i1];      // Copy from i2=n2-2 to i2=0
            }
        }

        // Copy faces in the n3 dimension (z-faces: i3=0 and i3=n3-1)
        // These copies are independent for each (i2, i1) pair.
        // Assuming n1, n2, n3 >= 3, these write locations (i3=0/n3-1)
        // do not overlap with reads/writes of the first two loops.
#pragma omp for // Parallelize over i2. Implicit barrier here syncs all three loops.
        for (i2 = 0; i2 < n2; i2++) {
            for (i1 = 0; i1 < n1; i1++) {
                u[n3 - 1][i2][i1] = u[1][i2][i1];    // Copy from i3=1 to i3=n3-1
                u[0][i2][i1] = u[n3 - 2][i2][i1];      // Copy from i3=n3-2 to i3=0
            }
        }

        // Implicit barrier at the end of the parallel region ensures all communication is complete.
    }
}