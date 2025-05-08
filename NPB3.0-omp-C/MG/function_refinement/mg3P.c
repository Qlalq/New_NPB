// Assuming necessary includes and definitions like lt, lb, m1, m2, m3, rprj3, zero3, psinv, interp, resid are available externally.

static void mg3P(double ****u, double ***v, double ****r, double a[4],
		 double c[4], int n1, int n2, int n3, int k) { // k is the level index parameter, which is modified internally
    int j;

    // --------------------------------------------------------------------
    // Loop 1: Restrict the residual (down the V-cycle)
    // This loop iterates over grid levels 'k' from lt (finest) down to lb+1.
    // This loop is inherently sequential across levels due to a loop-carried
    // data dependence: the computation of r at level k-1 depends on r at
    // level k. Specifically, r[k-1] is an output of the rprj3 function
    // which uses r[k] as an input.
    // Efficient OpenMP parallelism cannot be applied directly to this loop
    // over 'k'. Parallelism should be applied *within* the rprj3 function,
    // typically over the grid points (i, j, l) for a fixed level k.
    // --------------------------------------------------------------------
    for (k = lt; k >= lb+1; k--) { // k serves as the current level index
	j = k - 1; // j is the level index for the coarser result (one level below k)
	// Call rprj3 to restrict residual 'r' from level k to level j (k-1).
	// The signature implies: rprj3( output_r_level_j, ..., input_r_level_k, ... )
	rprj3(r[j], m1[k], m2[k], m3[k], // Output: r at level j
	      r[k], m1[j], m2[j], m3[j], k); // Input: r at level k; k is the current level
    }

    // --------------------------------------------------------------------
    // Block 2: Solve the equation on the coarsest grid (level lb)
    // These operations are performed solely at the coarsest level 'lb'.
    // They are sequential steps at this level (zero then solve).
    // Efficient OpenMP parallelism should be applied *within* the zero3
    // and psinv functions, typically over the grid points (i, j, l)
    // for level lb.
    // --------------------------------------------------------------------
    k = lb; // Set k to the coarsest level index 'lb'
    // Call zero3 to initialize the correction 'u' at level lb to zero.
    zero3(u[k], m1[k], m2[k], m3[k]);
    // Call psinv to perform the pseudo-inverse (coarse grid solve) at level lb.
    // This computes the approximate error 'u' from the residual 'r' at level lb.
    // Signature implies: psinv( input_r_level_k, output_u_level_k, ... )
    psinv(r[k], u[k], m1[k], m2[k], m3[k], c, k); // k is the current level lb

    // --------------------------------------------------------------------
    // Loop 3: Interpolate the error (up the V-cycle) and apply smoother
    // This loop iterates over grid levels 'k' from lb+1 up to lt-1.
    // This loop is inherently sequential across levels due to a loop-carried
    // data dependence: the computation of u at level k depends on u at
    // level k-1. Specifically, u[k] is an output of the interp function
    // which uses u[k-1] as an input.
    // Efficient OpenMP parallelism cannot be applied directly to this loop
    // over 'k'. Parallelism should be applied *within* the zero3, interp,
    // resid, and psinv functions, typically over the grid points (i, j, l)
    // for a fixed level k.
    // --------------------------------------------------------------------
    for (k = lb+1; k <= lt-1; k++) { // k serves as the current level index
	j = k - 1; // j is the level index for the coarser input (one level below k)

	// Call zero3 to re-initialize the correction 'u' at level k before interpolation.
	zero3(u[k], m1[k], m2[k], m3[k]);

	// Call interp to interpolate the error 'u' from level j (k-1) to level k,
	// and add it to the existing 'u' at level k.
	// Signature implies: interp( input_u_level_j, ..., output_u_level_k, ... )
	interp(u[j], m1[j], m2[j], m3[j], // Input: u at level j
	       u[k], m1[k], m2[k], m3[k], k); // Output: u at level k; k is the current level

	// Call resid to compute the residual 'r' using the updated 'u' at level k.
	// Signature implies: resid( input_u_level_k, input_r_level_k_or_v, output_r_level_k, ...)
	resid(u[k], r[k], r[k], m1[k], m2[k], m3[k], a, k); // k is the current level

	// Call psinv to apply the smoother (pseudo-inverse solve) at level k
	// based on the residual 'r' at level k. This updates 'u' at level k.
	// Signature implies: psinv( input_r_level_k, output_u_level_k, ... )
	psinv(r[k], u[k], m1[k], m2[k], m3[k], c, k); // k is the current level
    }

    // --------------------------------------------------------------------
    // Block 4: Final steps at the finest grid (level lt)
    // These operations are performed solely at the finest level 'lt'.
    // They are sequential steps at this level (interpolate, calc residual, solve).
    // The initial interp call depends on u[lt-1] computed in the last
    // iteration of Loop 3.
    // Efficient OpenMP parallelism should be applied *within* the interp,
    // resid, and psinv functions, typically over the grid points (i, j, l)
    // for level lt.
    // --------------------------------------------------------------------
    j = lt - 1; // j is the level index for the coarser input (lt-1)
    k = lt; // Set k to the finest level index 'lt'

    // Call interp to interpolate the error 'u' from level j (lt-1) to level k (lt),
    // and add it to the existing 'u' at level k.
    // Signature implies: interp( input_u_level_j, ..., output_u_level_k, ... )
    interp(u[j], m1[j], m2[j], m3[j], u[lt], n1, n2, n3, k); // k is the current level lt; use grid sizes n1, n2, n3

    // Call resid to compute the final residual 'r' using the updated 'u' at level lt
    // and the right-hand side 'v' at level lt.
    // Signature implies: resid( input_u_level_k, input_v, output_r_level_k, ...)
    resid(u[lt], v, r[lt], n1, n2, n3, a, k); // k is the current level lt; use grid sizes n1, n2, n3

    // Call psinv to apply the final smoother (pseudo-inverse solve) at level lt
    // based on the final residual 'r' at level lt. This updates 'u' at level lt.
    // Signature implies: psinv( input_r_level_k, output_u_level_k, ... )
    psinv(r[lt], u[lt], n1, n2, n3, c, k); // k is the current level lt; use grid sizes n1, n2, n3
}