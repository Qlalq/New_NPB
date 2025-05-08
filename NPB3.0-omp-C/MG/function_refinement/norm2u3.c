#include <stdio.h> // Retain original includes
#include <math.h>  // Required for fabs and sqrt

// Function refined for better structure and preparation for OpenMP
// - Separates the calculation of sum of squares and maximum absolute value
// - Uses the correct number of elements (interior points) for the norm calculation
static void norm2u3_refined(double ***r, int n1, int n2, int n3,
                            double *rnm2, double *rnmu, int nx, int ny, int nz) {

    // Variables to store the intermediate results of the two aggregations
    double sum_sq = 0.0;
    double max_abs = 0.0;

    int i3, i2, i1; // Loop indices

    // Determine the number of elements being included in the calculations (interior points)
    // Using long long to prevent potential overflow if dimensions are large
    long long num_interior_elements = 0;
    if (n1 > 2 && n2 > 2 && n3 > 2) {
        num_interior_elements = (long long)(n1-2) * (long long)(n2-2) * (long long)(n3-2);
    } else {
        // If dimensions are too small, there are no interior points
        num_interior_elements = 0;
    }

    // --- Decouple and Structure: Separate the two reduction operations ---
    // This loop calculates the sum of squares of the interior elements.
    // It is prepared for a parallel loop with a sum reduction clause.
    // #pragma omp parallel for collapse(3) reduction(+:sum_sq)
    for (i3 = 1; i3 < n3-1; i3++) {
        for (i2 = 1; i2 < n2-1; i2++) {
            for (i1 = 1; i1 < n1-1; i1++) {
                double element_val = r[i3][i2][i1];
                sum_sq += element_val * element_val; // Accumulation
            }
        }
    }

    // This loop calculates the maximum absolute value of the interior elements.
    // It is prepared for a parallel loop with a maximum reduction clause.
    // #pragma omp parallel for collapse(3) reduction(max:max_abs)
    for (i3 = 1; i3 < n3-1; i3++) {
        for (i2 = 1; i2 < n2-1; i2++) {
            for (i1 = 1; i1 < n1-1; i1++) {
                double a = fabs(r[i3][i2][i1]);
                if (a > max_abs) max_abs = a; // Maximum
            }
        }
    }
    // --- End Decouple and Structure ---


    // --- Final calculations and assignment ---
    // Assign the maximum absolute value
    *rnmu = max_abs;

    // Calculate and assign the L2 norm (RMS)
    // Handle the case where there are no interior elements to avoid division by zero
    if (num_interior_elements > 0) {
        *rnm2 = sqrt(sum_sq / (double)num_interior_elements);
    } else {
        // If no interior elements, the norm is 0
        *rnm2 = 0.0;
    }
    // --- End Final calculations ---
}