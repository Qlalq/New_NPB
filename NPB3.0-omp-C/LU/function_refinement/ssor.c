static void ssor(void) {
  int i, j, k, m;
  int istep;
  double  tmp;
  // tv is likely a temporary array used within buts.
  // Its size is related to the spatial grid dimensions ISIZ1, ISIZ2.
  double  delunm[5], tv[ISIZ1][ISIZ2][5];
  tmp = 1.0 / ( omega * ( 2.0 - omega ) ) ;

  // Initialization Phase:
  // This block initializes arrays a, b, c, d to zero.
  // The nested loops are embarrassingly parallel. Each iteration
  // writes to a unique element based on the loop indices (i, j, k, m),
  // making them independent. This is a good candidate for parallelization
  // over any of the outer dimensions (i, j, or k).
  {
    for (i = 0; i < ISIZ1; i++) {
      for (j = 0; j < ISIZ2; j++) {
	for (k = 0; k < 5; k++) { // Note: k loop here is 0-4, array dimension, not spatial k
	  for (m = 0; m < 5; m++) { // Note: m loop here is 0-4, array dimension, not spatial m
	    a[i][j][k][m] = 0.0;
	    b[i][j][k][m] = 0.0;
	    c[i][j][k][m] = 0.0;
	    d[i][j][k][m] = 0.0;
	  }
	}
      }
    }
  }

  // rhs computation:
  // This function likely computes the right-hand side vector.
  // It typically involves grid point operations. Parallelization
  // opportunities are expected within this function, likely over
  // the spatial dimensions (i, j, k).
  rhs();

  // l2norm computation:
  // This function computes the L2 norm, which involves summing up
  // values across the domain and then taking a square root. The summation
  // is a reduction operation, which can be efficiently parallelized using
  // OpenMP reduction clauses.
  l2norm( nx0, ny0, nz0,
	  ist, iend, jst, jend,
	  rsd, rsdnm );

  timer_clear(1);
  timer_start(1);

  // Main time stepping loop (SSOR iteration):
  // The SSOR method is iterative, meaning the computation in step 'istep'
  // depends on the results from step 'istep-1'. Therefore, the outer
  // 'istep' loop is inherently sequential and cannot be parallelized directly
  // without changing the algorithm (e.g., to a different solver).
  // Parallelization efforts should focus on the computations *within* this loop.
  for (istep = 1; istep <= itmax; istep++) {

    if (istep%20  ==  0 || istep  ==  itmax || istep  ==  1) {
      printf(" Time step %4d\n", istep); // Print statement is sequential
    }

    // Operations within a single time step:
    // This block contains the computations that are performed for each
    // time step. These computations are the primary targets for parallelization.
    {
      // Scaling rsd:
      // This nested loop scales the values in the rsd array by dt.
      // Similar to the initialization loops, each iteration (i, j, k, m)
      // updates a unique element of rsd based only on its previous value.
      // There are no loop-carried dependencies across i, j, k, or m within
      // this block. This is embarrassingly parallel and can be parallelized
      // over any of the spatial dimensions (i, j, or k).
      for (i = ist; i <= iend; i++) {
	for (j = jst; j <= jend; j++) {
	  for (k = 1; k <= nz - 2; k++) {
	    for (m = 0; m < 5; m++) {
	      rsd[i][j][k][m] = dt * rsd[i][j][k][m];
	    }
	  }
	}
      }

      // Forward sweep (BLTS - Block Lower Triangular Solve):
      // This loop performs a forward substitution along the k-dimension.
      // The computation for k depends on the results from previous k values.
      // Therefore, the loop over k is typically sequential due to data dependencies.
      // However, the operations performed *within* jacld and blts for a fixed k
      // often involve computations over the i and j dimensions (a 2D plane).
      // The primary parallelization opportunity here is expected to be found
      // by parallelizing the operations *within* the jacld and blts functions,
      // likely over i and j for each k slice.
      for (k = 1; k <= nz - 2; k++) {
	jacld(k); // Likely computes matrix coefficients for plane k
	blts(nx, ny, nz, k, omega, rsd, a, b, c, d, ist, iend, jst, jend, nx0, ny0 ); // Solves for plane k
      }

      // Backward sweep (BUTS - Block Upper Triangular Solve):
      // This loop performs a backward substitution along the k-dimension.
      // The computation for k depends on the results from subsequent k values.
      // Similar to the forward sweep, the loop over k is typically sequential.
      // Parallelization opportunities are expected *within* the jacu and buts
      // functions, likely over i and j for each k slice.
      for (k = nz - 2; k >= 1; k--) {
	jacu(k); // Likely computes matrix coefficients for plane k
	buts(nx, ny, nz, k, omega, rsd, tv, d, a, b, c, ist, iend, jst, jend, nx0, ny0 ); // Solves for plane k
      }

      // Updating solution u:
      // This nested loop updates the solution array u using the computed residual rsd.
      // Each iteration (i, j, k, m) updates a unique element of u based on its
      // previous value and the corresponding element of rsd. There are no
      // loop-carried dependencies across i, j, k, or m within this block.
      // This is embarrassingly parallel and can be parallelized over any
      // of the spatial dimensions (i, j, or k).
      for (i = ist; i <= iend; i++) {
	for (j = jst; j <= jend; j++) {
	  for (k = 1; k <= nz-2; k++) {
	    for (m = 0; m < 5; m++) {
	      u[i][j][k][m] = u[i][j][k][m]
		+ tmp * rsd[i][j][k][m];
	    }
	  }
	}
      }
    } // End of operations within a time step

    // Conditional l2norm computation:
    // Another L2 norm calculation. Involves reduction, parallelizable.
    if ( istep % inorm  ==  0 ) {
      l2norm( nx0, ny0, nz0,
	      ist, iend, jst, jend,
	      rsd, delunm );
    }

    // rhs computation: Potentially parallelizable internally over spatial dimensions.
    rhs();

    // Conditional l2norm computation: Involves reduction, parallelizable.
    if ( ( istep % inorm  ==  0 ) ||
	 ( istep  ==  itmax ) ) {
      l2norm( nx0, ny0, nz0,
	      ist, iend, jst, jend,
	      rsd, rsdnm );
    }

    // Convergence check: This check compares elements of the rsdnm array
    // to tolerance values. This comparison logic is sequential, but the
    // rsdnm array itself is the result of a parallelizable reduction.
    if ( ( rsdnm[0] < tolrsd[0] ) &&
	 ( rsdnm[1] < tolrsd[1] ) &&
	 ( rsdnm[2] < tolrsd[2] ) &&
	 ( rsdnm[3] < tolrsd[3] ) &&
	 ( rsdnm[4] < tolrsd[4] ) ) {
	exit(1); // Program exit is sequential
    }
  } // End of time stepping loop

  timer_stop(1); // Sequential timer stop
  maxtime= timer_read(1); // Sequential read
}