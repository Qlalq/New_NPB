static void resid_refined( double ***u, double ***v, double ***r,
		   int n1, int n2, int n3, double a[4], int k ) {
    int i3, i2, i1;

    // The loops iterate over a 3D grid, calculating residuals.
    // The computation for each point r[i3][i2][i1] depends on a stencil
    // of points in u and v. Critically, the calculation of one r[i3][i2][i1]
    // does NOT depend on the value of r at any other point.
    // The outer loops (i3 and i2) define independent 'planes' of computation.
    // This structure is highly suitable for parallelization over i3 or i2.

    for (i3 = 1; i3 < n3-1; i3++) {
	for (i2 = 1; i2 < n2-1; i2++) {

            // The temporary arrays u1 and u2 store intermediate results
            // for the current (i3, i2) plane.
            // By declaring them inside the i2 loop, they become local to each
            // (i3, i2) iteration. When parallelizing the i3 or i2 loop with OpenMP,
            // these arrays can be easily marked as private for each thread
            // or for each chunk of iterations assigned to a thread.
            // This avoids false sharing and minimizes synchronization overhead
            // related to these temporaries.
            // Assuming M is a constant or macro defined elsewhere, typically >= n1.
            double u1[M], u2[M];

            // First inner loop: Pre-calculate stencil parts for the current (i3, i2) plane.
            // This loop has no loop-carried dependencies. All iterations are independent
            // and write to unique elements of the temporary arrays u1 and u2.
            for (i1 = 0; i1 < n1; i1++) {
		u1[i1] = u[i3][i2-1][i1] + u[i3][i2+1][i1]
		       + u[i3-1][i2][i1] + u[i3+1][i2][i1];
		u2[i1] = u[i3-1][i2-1][i1] + u[i3-1][i2+1][i1]
		       + u[i3+1][i2-1][i1] + u[i3+1][i2+1][i1];
	    }

            // Second inner loop: Calculate and write the residual r for the current (i3, i2) plane.
            // This loop writes to r[i3][i2][i1]. For a fixed (i3, i2), each i1 writes to a unique
            // element of r. There are no loop-carried dependencies on r within this loop.
            // It reads from the u1 and u2 arrays computed in the previous loop for the same (i3, i2).
	    for (i1 = 1; i1 < n1-1; i1++) {
		r[i3][i2][i1] = v[i3][i2][i1]
		    - a[0] * u[i3][i2][i1]
		    - a[2] * ( u2[i1] + u1[i1-1] + u1[i1+1] ) // Reads from u1 and u2
		    - a[3] * ( u2[i1-1] + u2[i1+1] );   // Reads from u1 and u2
	    }
	}
    }

    // Subsequent calls are outside the main computation loops.
    // These might involve communication or debugging/reporting.
    // Their parallelization or synchronization needs are separate from the
    // stencil computation loops refined above.
    comm3(r,n1,n2,n3,k);
    if (debug_vec[0] >= 1 ) {
	rep_nrm(r,n1,n2,n3,"   resid",k);
    }
    if ( debug_vec[2] >= k ) {
	showall(r,n1,n2,n3);
    }
}