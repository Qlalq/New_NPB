#include <stdio.h> // Required for printf

// Define min if it's not available via other includes
// This ensures the code is self-contained for the snippet.
#ifndef min
#define min(a,b) ((a) < (b) ? (a) : (b))
#endif

// Refined function structure for potential parallelization analysis.
// Note: While the data accesses are independent, the use of printf for
// formatted sequential output inherently introduces dependencies on the
// output stream order, severely limiting the practical efficiency of
// parallelizing these specific loops for this purpose.
// This refinement primarily focuses on clearer variable scope and explicit data access.
static void showall_refined(double ***z, int n1, int n2, int n3) {
    // Determine the actual limits for printing based on n1, n2, n3 and maximum display size.
    // These limits are calculated once before the loops, avoiding recalculation.
    int m1 = min(n1, 18);
    int m2 = min(n2, 14);
    int m3 = min(n3, 18);

    printf("\n"); // Sequential I/O operation before the main output block

    // The loops iterate through the data blocks and elements to print.
    // The order of iteration (i3 -> i1 -> i2) is crucial for producing
    // the required structured output format (blocks -> rows -> elements within rows).
    // Each iteration of the innermost loop reads an independent data element z[i3][i2][i1].
    // However, the primary operation is printf, which writes to a shared resource (stdout).
    // To maintain the specific output order and formatting, these printf calls
    // must execute mostly sequentially or require complex synchronization,
    // making direct parallelization for performance gain impractical for this I/O bound task.

    // The loop counters (i3, i1, i2) are declared within the loop scope.
    // This is a standard C practice that limits variable scope and is beneficial
    // for parallelization as loop counters are typically privatized by each thread.

    // Iterate through the outer dimension (blocks) - Potential parallelization target
    // if inner work was independent and computationally significant, ignoring I/O.
    for (int i3 = 0; i3 < m3; ++i3) {

        // Iterate through the middle dimension (rows within blocks)
        for (int i1 = 0; i1 < m1; ++i1) {

            // Iterate through the inner dimension (elements within rows)
            // Each access to z[i3][i2][i1] is independent of other iterations.
            // This loop must execute sequentially or with strict synchronization
            // to ensure elements within a row are printed in the correct order.
            for (int i2 = 0; i2 < m2; ++i2) {
                // Explicitly separate the data access operation (reading from z)
                // from the side-effect operation (printing).
                double value_to_print = z[i3][i2][i1];

                // Perform the sequential I/O operation for the current element.
                // This operation is the main bottleneck for parallelization here.
                printf("%6.3f", value_to_print);
            }
            // Sequential I/O operation for the newline after each row.
            printf("\n");
        }
        // Sequential I/O operation for the separator between blocks.
        printf(" - - - - - - - \n");
    }

    printf("\n"); // Sequential I/O operation after the main output block
}