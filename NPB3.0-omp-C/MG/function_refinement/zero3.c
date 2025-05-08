static void zero3(double ***z, int n1, int n2, int n3) {
    int i1, i2, i3;
    // This structure is already well-suited for OpenMP parallelization.
    // 1. Data Dependencies: Each element z[i3][i2][i1] is written independently.
    //    There are no loop-carried dependencies.
    // 2. Load Imbalance: The work inside the innermost loop (assignment) is constant.
    //    Parallelizing any of the outer loops (i3 or i2) provides good granularity
    //    and ensures even distribution of work among threads.
    // 3. Synchronization Overhead: Since each thread can operate on a distinct portion
    //    of the array without accessing memory being written by another thread,
    //    explicit synchronization for the write operation itself is not required
    //    when using constructs like #pragma omp parallel for.
    for (i3 = 0; i3 < n3; i3++) {
	for (i2 = 0; i2 < n2; i2++) {
            for (i1 = 0; i1 < n1; i1++) {
		z[i3][i2][i1] = 0.0;
	    }
	}
    }
}