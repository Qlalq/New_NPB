static void ssor(void) {
  int i, j, k, m;
  int istep;
  double  tmp;
  double  delunm[5], tv[ISIZ1][ISIZ2][5];
  tmp = 1.0 / ( omega * ( 2.0 - omega ) ) ;
{
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
  rhs();
  l2norm( nx0, ny0, nz0,
	  ist, iend, jst, jend,
	  rsd, rsdnm );
  timer_clear(1);
  timer_start(1);
  for (istep = 1; istep <= itmax; istep++) {
    if (istep%20  ==  0 || istep  ==  itmax || istep  ==  1) {
      printf(" Time step %4d\n", istep);
    }
{  
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
      jacld(k);
      blts(nx, ny, nz, k,
	   omega,
	   rsd,
	   a, b, c, d,
	   ist, iend, jst, jend, 
	   nx0, ny0 );
    }
    for (k = nz - 2; k >= 1; k--) {
      jacu(k);
      buts(nx, ny, nz, k,
	   omega,
	   rsd, tv,
	   d, a, b, c,
	   ist, iend, jst, jend,
	   nx0, ny0 );
    }
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
} 
    if ( istep % inorm  ==  0 ) {
      l2norm( nx0, ny0, nz0,
	      ist, iend, jst, jend,
	      rsd, delunm );
    }
    rhs();
    if ( ( istep % inorm  ==  0 ) ||
	 ( istep  ==  itmax ) ) {
      l2norm( nx0, ny0, nz0,
	      ist, iend, jst, jend,
	      rsd, rsdnm );
    }
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