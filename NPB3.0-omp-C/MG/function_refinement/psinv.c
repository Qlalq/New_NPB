// Note: The constant M must be defined elsewhere (e.g., #define M ...)
// and its value must be sufficient (M >= n1) for the array declarations.
// Necessary includes (like stdio.h, math.h depending on comm3, rep_nrm, showall, etc.)
// and declarations (like debug_vec, comm3, rep_nrm, showall) are assumed to be
// available from the surrounding code base.

static void psinv( double ***r, double ***u, int n1, int n2, int n3,
		   double c[4], int k) {
    int i3, i2, i1;

    // The outermost loops (i3 and i2) iterate over the interior grid points.
    // These loops are the primary candidates for parallelization as their
    // iterations operate on largely independent parts of the grid (different
    // i3-planes or i2-planes).
    for (i3 = 1; i3 < n3-1; i3++) {
	for (i2 = 1; i2 < n2-1; i2++) {
            // Declare temporary arrays r1 and r2 here. This scopes them to be
            // local to each iteration of the i2 loop. In OpenMP, variables
            // declared within the loop body of a parallelized loop (like i2
            // or i3 if i2 is also parallelized/nested) are typically private
            // by default. This structure explicitly shows that r1 and r2 are
            // temporary scratch space needed per (i3, i2) work unit.
            // This addresses the privatization requirement for efficient parallelization
            // by avoiding false sharing or race conditions on these temporary buffers.
            double r1[M], r2[M];

            // First inner loop: Compute intermediate stencil values r1 and r2.
            // This loop iterates n1 times. There are no loop-carried dependencies
            // within this loop; each element r1[i1] and r2[i1] is computed based
            // on values from 'r' at fixed offsets relative to the current (i3, i2, i1).
            for (i1 = 0; i1 < n1; i1++) {
		r1[i1] = r[i3][i2-1][i1] + r[i3][i2+1][i1]
		    + r[i3-1][i2][i1] + r[i3+1][i2][i1];
		r2[i1] = r[i3-1][i2-1][i1] + r[i3-1][i2+1][i1]
		    + r[i3+1][i2-1][i1] + r[i3+1][i2+1][i1];
	    }

            // Second inner loop: Update the 'u' array using a stencil computation.
            // This computation uses values from 'r' and the intermediate values
            // r1 and r2 computed in the previous loop for the same (i3, i2) pair.
            // This loop iterates n1-2 times. There are no loop-carried dependencies
            // on the write destination u[i3][i2][i1]; each iteration writes to a
            // unique element of the 'u' array within this (i3, i2) plane.
            // Dependencies on r1[i1-1] and r1[i1+1] are read-only access to data
            // that was fully computed and stored in r1 during the previous loop.
            for (i1 = 1; i1 < n1-1; i1++) {
		u[i3][i2][i1] = u[i3][i2][i1]
		    + c[0] * r[i3][i2][i1]
		    + c[1] * ( r[i3][i2][i1-1] + r[i3][i2][i1+1]
			       + r1[i1] )
		    + c[2] * ( r2[i1] + r1[i1-1] + r1[i1+1] );
	    }
            // The sequential execution of the two i1 loops within each (i3, i2)
            // iteration ensures that r1 and r2 are fully populated before being used.
            // The independent nature of updates to u[i3][i2][i1] for different
            // (i3, i2) pairs allows for parallelization of the outer loops.
	} // r1 and r2 go out of scope here, their storage can be reused by other threads/iterations
    }

    // Communication step. This function likely performs boundary exchanges
    // or other collective operations on the 'u' array. It must be called
    // after the core updates to the interior points are complete across all
    // parallel tasks. It typically acts as a synchronization point.
    comm3(u,n1,n2,n3,k);

    // Debugging and reporting steps. These functions read the results from 'u'.
    // They should be safe to call after the communication step ensures 'u' is
    // consistent, but could potentially be performance bottlenecks depending
    // on their implementation (e.g., I/O, collective norms).
    if (debug_vec[0] >= 1 ) {
	rep_nrm(u,n1,n2,n3,"   psinv",k);
    }
    if ( debug_vec[3] >= k ) {
	showall(u,n1,n2,n3);
    }
}