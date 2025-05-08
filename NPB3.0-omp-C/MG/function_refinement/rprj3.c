static void rprj3( double ***r, int m1k, int m2k, int m3k,
		   double ***s, int m1j, int m2j, int m3j, int k ) {
    int j3, j2, j1;
    int i3, i2, i1;
    int d1, d2, d3; // d variables can be calculated once as they depend only on input parameters

    // Calculate d1, d2, d3 outside the main loops
    if (m1k == 3) {
        d1 = 2;
    } else {
        d1 = 1;
    }
    if (m2k == 3) {
        d2 = 2;
    } else {
        d2 = 1;
    }
    if (m3k == 3) {
        d3 = 2;
    } else {
        d3 = 1;
    }

    // The outer j3 loop. This loop's iterations are independent.
    // Could be parallelized, possibly collapsed with the j2 loop.
    for (j3 = 1; j3 < m3j-1; j3++) {
	i3 = 2*j3-d3;

        // The j2 loop. This loop's iterations are independent from each other,
        // provided the temporary arrays x1 and y1 are handled correctly.
        // This is a primary candidate for OpenMP parallelization.
        for (j2 = 1; j2 < m2j-1; j2++) {
            i2 = 2*j2-d2;

            // Declare temporary arrays x1, y1 and scalars x2, y2 inside
            // the j2 loop scope. This ensures that each iteration of j2
            // (and thus each thread processing a j2 chunk when parallelized)
            // gets its own private copy of these temporaries.
            // This eliminates data dependencies between j2 iterations
            // that were previously caused by x1, y1, x2, y2 being shared.
            // M is assumed to be a predefined constant large enough
            // for indices 2*j1-d1 up to 2*(m1j-1)-d1 + 2.
            double x1[M], y1[M];
            double x2, y2;

            // First inner loop j1: Calculate x1 and y1 values
            // These calculations are independent for each j1 iteration within this j2 block.
            // This loop can run sequentially within the j2 task, or could be parallelized
            // internally if the cost is high and M is large.
            for (j1 = 1; j1 < m1j; j1++) { // Original loop goes up to m1j-1 (exclusive) -> [1, m1j-1]
                i1 = 2*j1-d1;
                x1[i1] = r[i3+1][i2][i1] + r[i3+1][i2+2][i1]
                    + r[i3][i2+1][i1] + r[i3+2][i2+1][i1];
                y1[i1] = r[i3][i2][i1] + r[i3+2][i2][i1]
                    + r[i3][i2+2][i1] + r[i3+2][i2+2][i1];
            }

            // Second inner loop j1: Calculate s using r and the temporary x1, y1 values
            // These calculations are independent for each j1 iteration within this j2 block.
            // This loop depends on the results of the first j1 loop for this j2 block.
            // This loop can run sequentially within the j2 task, or could be parallelized
            // internally.
            for (j1 = 1; j1 < m1j-1; j1++) { // Original loop goes up to m1j-1 (exclusive) -> [1, m1j-2]
                i1 = 2*j1-d1;
                y2 = r[i3][i2][i1+1] + r[i3+2][i2][i1+1]
                    + r[i3][i2+2][i1+1] + r[i3+2][i2+2][i1+1];
                x2 = r[i3+1][i2][i1+1] + r[i3+1][i2+2][i1+1]
                    + r[i3][i2+1][i1+1] + r[i3+2][i2+1][i1+1];
                s[j3][j2][j1] =
                    0.5 * r[i3+1][i2+1][i1+1]
                    + 0.25 * ( r[i3+1][i2+1][i1] + r[i3+1][i2+1][i1+2] + x2)
                    + 0.125 * ( x1[i1] + x1[i1+2] + y2) // Reads from x1, y1
                    + 0.0625 * ( y1[i1] + y1[i1+2] );   // Reads from x1, y1
            }
	}
    }
    // Post-processing steps - usually remain sequential or require separate parallelization/synchronization
    comm3(s,m1j,m2j,m3j,k-1);
    if (debug_vec[0] >= 1 ) {
	rep_nrm(s,m1j,m2j,m3j,"   rprj3",k-1);
    }
    if (debug_vec[4] >= k ) {
	showall(s,m1j,m2j,m3j);
    }
}