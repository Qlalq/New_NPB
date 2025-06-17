static void ssor(void) {

/*--------------------------------------------------------------------
c   to perform pseudo-time stepping SSOR iterations
c   for five nonlinear pde s.
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j, k, m;
  int istep;
  double  tmp;
  double  delunm[5], tv[ISIZ1][ISIZ2][5];

/*--------------------------------------------------------------------
c   begin pseudo-time stepping iterations
--------------------------------------------------------------------*/
  tmp = 1.0 / ( omega * ( 2.0 - omega ) ) ;

/*--------------------------------------------------------------------
c   initialize a,b,c,d to zero (guarantees that page tables have been
c   formed, if applicable on given architecture, before timestepping).
--------------------------------------------------------------------*/
#pragma omp parallel private(i,j,k,m)
{
#pragma omp for    
  for (i = 0; i < ISIZ1; i++) {
    for (j = 0; j < ISIZ2; j++) {
      for (k = 0; k < 5; k++) {
	for (m = 0; m < 5; m++) {
	  a[i][j][k][m] = 0.0;
	  b[i][j][k][m] = 0.0;
	  c[i][j][k][m] = 0.0;
	  d[i][j][k][m] = 0.0;
	}
      }
    }
  }
}
/*--------------------------------------------------------------------
c   compute the steady-state residuals
--------------------------------------------------------------------*/
  rhs();
 
/*--------------------------------------------------------------------
c   compute the L2 norms of newton iteration residuals
--------------------------------------------------------------------*/
  l2norm( nx0, ny0, nz0,
	  ist, iend, jst, jend,
	  rsd, rsdnm );


  timer_clear(1);
  timer_start(1);
 
/*--------------------------------------------------------------------
c   the timestep loop
--------------------------------------------------------------------*/


  for (istep = 1; istep <= itmax; istep++) {

    if (istep%20  ==  0 || istep  ==  itmax || istep  ==  1) {
#pragma omp master	
      printf(" Time step %4d\n", istep);
    }

#pragma omp parallel private(istep,i,j,k,m)
{  
 
/*--------------------------------------------------------------------
c   perform SSOR iteration
--------------------------------------------------------------------*/
#pragma omp for    
    for (i = ist; i <= iend; i++) {
      for (j = jst; j <= jend; j++) {
	for (k = 1; k <= nz - 2; k++) {
	  for (m = 0; m < 5; m++) {
	    rsd[i][j][k][m] = dt * rsd[i][j][k][m];
	  }
	}
      }
    }

    for (k = 1; k <= nz - 2; k++) {
/*--------------------------------------------------------------------
c   form the lower triangular part of the jacobian matrix
--------------------------------------------------------------------*/
      jacld(k);
 
/*--------------------------------------------------------------------
c   perform the lower triangular solution
--------------------------------------------------------------------*/
      blts(nx, ny, nz, k,
	   omega,
	   rsd,
	   a, b, c, d,
	   ist, iend, jst, jend, 
	   nx0, ny0 );
    }
    
#pragma omp barrier

    for (k = nz - 2; k >= 1; k--) {
/*--------------------------------------------------------------------
c   form the strictly upper triangular part of the jacobian matrix
--------------------------------------------------------------------*/
      jacu(k);

/*--------------------------------------------------------------------
c   perform the upper triangular solution
--------------------------------------------------------------------*/
      buts(nx, ny, nz, k,
	   omega,
	   rsd, tv,
	   d, a, b, c,
	   ist, iend, jst, jend,
	   nx0, ny0 );
    }
#pragma omp barrier 
 
/*--------------------------------------------------------------------
c   update the variables
--------------------------------------------------------------------*/

#pragma omp for
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
} /* end parallel */
/*--------------------------------------------------------------------
c   compute the max-norms of newton iteration corrections
--------------------------------------------------------------------*/
    if ( istep % inorm  ==  0 ) {
      l2norm( nx0, ny0, nz0,
	      ist, iend, jst, jend,
	      rsd, delunm );
    }
 
/*--------------------------------------------------------------------
c   compute the steady-state residuals
--------------------------------------------------------------------*/
    rhs();
 
/*--------------------------------------------------------------------
c   compute the max-norms of newton iteration residuals
--------------------------------------------------------------------*/
    if ( ( istep % inorm  ==  0 ) ||
	 ( istep  ==  itmax ) ) {
      l2norm( nx0, ny0, nz0,
	      ist, iend, jst, jend,
	      rsd, rsdnm );
    }

/*--------------------------------------------------------------------
c   check the newton-iteration residuals against the tolerance levels
--------------------------------------------------------------------*/
    if ( ( rsdnm[0] < tolrsd[0] ) &&
	 ( rsdnm[1] < tolrsd[1] ) &&
	 ( rsdnm[2] < tolrsd[2] ) &&
	 ( rsdnm[3] < tolrsd[3] ) &&
	 ( rsdnm[4] < tolrsd[4] ) ) {
	exit(1);
    }
  }

 
  timer_stop(1);
  maxtime= timer_read(1);
 
}