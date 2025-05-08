static void setup(int *n1, int *n2, int *n3, int lt) {
    int k;

    // Section 1: Compute grid sizes for coarser levels.
    // This loop has a loop-carried dependency (nx[k] depends on nx[k+1]),
    // making it inherently sequential in this form.
    // OpenMP parallelization is not straightforward for this loop alone
    // without a different algorithm or approach.
    for ( k = lt-1; k >= 1; k--) {
        nx[k] = nx[k+1]/2;
        ny[k] = ny[k+1]/2;
        nz[k] = nz[k+1]/2;
    }

    // Section 2: Compute dimension sizes m1, m2, m3 based on grid sizes.
    // The calculations for each index 'k' are independent of each other.
    // This loop is well-structured for OpenMP parallelization.
    // Each thread can compute a different range of 'k' values without
    // requiring synchronization within the loop body, as writes to m1, m2, m3
    // at index k do not conflict with other indices.
    for (k = 1; k <= lt; k++) {
        m1[k] = nx[k]+2;
        m2[k] = nz[k]+2; // Uses nz based on original code.
        m3[k] = ny[k]+2;
    }

    // Section 3: Set boundary conditions and output sizes for the finest level (lt).
    // These are scalar assignments depending on values computed in Section 1.
    // They are typically executed sequentially after the grid sizes are determined.
    // No parallelization benefit for these few assignments.
    is1 = 1;
    ie1 = nx[lt];
    *n1 = nx[lt]+2;

    is2 = 1;
    ie2 = ny[lt];
    *n2 = ny[lt]+2;

    is3 = 1;
    ie3 = nz[lt];
    *n3 = nz[lt]+2;

    // Section 4: Optional debug print statement.
    // This is a sequential output operation.
    // If executed within a parallel region, it would require synchronization
    // (e.g., #pragma omp critical or #pragma omp master) to avoid mixed output,
    // but is logically outside the main computational loops.
    if (debug_vec[1] >=  1 ) {
        printf(" in setup, \n");
        printf("  lt  nx  ny  nz  n1  n2  n3 is1 is2 is3 ie1 ie2 ie3\n");
        printf("%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d\n",
               lt,nx[lt],ny[lt],nz[lt],*n1,*n2,*n3,is1,is2,is3,ie1,ie2,ie3);
    }
}