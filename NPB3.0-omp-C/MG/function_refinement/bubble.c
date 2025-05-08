#define M /* Assuming M is defined elsewhere, keeping it as in original */

static void bubble_refined( double ten[M][2], int j1[M][2], int j2[M][2],
		    int j3[M][2], int m, int ind ) {
    double temp;
    int i, j_temp;

    /*
     * Refinement for OpenMP:
     * 1. Data Dependencies: The original loop iterates i from 0 to m-2, comparing/swapping (i, i+1).
     *    This creates a loop-carried dependency between i and i+1. To break this,
     *    we split the single pass into two sub-passes: one for odd-indexed pairs
     *    (1,2), (3,4), ... and one for even-indexed pairs (0,1), (2,3), ...
     *    Iterations within each sub-pass are independent.
     * 2. Load Imbalance: The original code's 'return' statement caused unpredictable
     *    early exit, leading to severe load imbalance if parallelized naively.
     *    Removing the 'return' ensures a fixed amount of work per function call
     *    (covering all adjacent pairs). Splitting into odd/even passes distributes
     *    this work across two loops, which can be parallelized for better balance.
     * 3. Synchronization Overhead: By making iterations within each odd/even loop
     *    independent, a simple '#pragma omp for' can be applied to each loop
     *    without needing fine-grained synchronization (like locks) on individual
     *    array elements, reducing overhead compared to attempting to parallelize
     *    the original dependent loop structure.
     */

    /*
     * Pass over odd-indexed pairs (1,2), (3,4), (5,6), ...
     * Loop i from 1 up to m-2, stepping by 2.
     * This loop can be parallelized using OpenMP.
     */
    for (i = 1; i < m - 1; i += 2) {
        if (ind == 1) { /* Descending sort for ten[][ind] */
            if ( ten[i][ind] > ten[i+1][ind] ) {
                /* Swap ten[i][ind] and ten[i+1][ind] */
                temp = ten[i+1][ind];
                ten[i+1][ind] = ten[i][ind];
                ten[i][ind] = temp;
                /* Swap corresponding elements in j1[][ind] */
                j_temp = j1[i+1][ind];
                j1[i+1][ind] = j1[i][ind];
                j1[i][ind] = j_temp;
                /* Swap corresponding elements in j2[][ind] */
                j_temp = j2[i+1][ind];
                j2[i+1][ind] = j2[i][ind];
                j2[i][ind] = j_temp;
                /* Swap corresponding elements in j3[][ind] */
                j_temp = j3[i+1][ind];
                j3[i+1][ind] = j3[i][ind];
                j3[i][ind] = j_temp;
            }
        } else { /* Ascending sort for ten[][ind] */
            if ( ten[i][ind] < ten[i+1][ind] ) {
                /* Swap ten[i][ind] and ten[i+1][ind] */
                temp = ten[i+1][ind];
                ten[i+1][ind] = ten[i][ind];
                ten[i][ind] = temp;
                /* Swap corresponding elements in j1[][ind] */
                j_temp = j1[i+1][ind];
                j1[i+1][ind] = j1[i][ind];
                j1[i][ind] = j_temp;
                /* Swap corresponding elements in j2[][ind] */
                j_temp = j2[i+1][ind];
                j2[i+1][ind] = j2[i][ind];
                j2[i][ind] = j_temp;
                /* Swap corresponding elements in j3[][ind] */
                j_temp = j3[i+1][ind];
                j3[i+1][ind] = j3[i][ind];
                j3[i][ind] = j_temp;
            }
        }
    }

    /*
     * Pass over even-indexed pairs (0,1), (2,3), (4,5), ...
     * Loop i from 0 up to m-2, stepping by 2.
     * This loop can be parallelized using OpenMP.
     */
    for (i = 0; i < m - 1; i += 2) {
        if (ind == 1) { /* Descending sort for ten[][ind] */
            if ( ten[i][ind] > ten[i+1][ind] ) {
                /* Swap ten[i][ind] and ten[i+1][ind] */
                temp = ten[i+1][ind];
                ten[i+1][ind] = ten[i][ind];
                ten[i][ind] = temp;
                /* Swap corresponding elements in j1[][ind] */
                j_temp = j1[i+1][ind];
                j1[i+1][ind] = j1[i][ind];
                j1[i][ind] = j_temp;
                /* Swap corresponding elements in j2[][ind] */
                j_temp = j2[i+1][ind];
                j2[i+1][ind] = j2[i][ind];
                j2[i][ind] = j_temp;
                /* Swap corresponding elements in j3[][ind] */
                j_temp = j3[i+1][ind];
                j3[i+1][ind] = j3[i][ind];
                j3[i][ind] = j_temp;
            }
        } else { /* Ascending sort for ten[][ind] */
            if ( ten[i][ind] < ten[i+1][ind] ) {
                /* Swap ten[i][ind] and ten[i+1][ind] */
                temp = ten[i+1][ind];
                ten[i+1][ind] = ten[i][ind];
                ten[i][ind] = temp;
                /* Swap corresponding elements in j1[][ind] */
                j_temp = j1[i+1][ind];
                j1[i+1][ind] = j1[i][ind];
                j1[i][ind] = j_temp;
                /* Swap corresponding elements in j2[][ind] */
                j_temp = j2[i+1][ind];
                j2[i+1][ind] = j2[i][ind];
                j2[i][ind] = j_temp;
                /* Swap corresponding elements in j3[][ind] */
                j_temp = j3[i+1][ind];
                j3[i+1][ind] = j3[i][ind];
                j3[i][ind] = j_temp;
            }
        }
    }
}