static void resid( double ***u, double ***v, double ***r,
		   int n1, int n2, int n3, double a[4], int k ) {

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c     resid computes the residual:  r = v - Au
c
c     This  implementation costs  15A + 4M per result, where
c     A and M denote the costs of Addition (or Subtraction) and 
c     Multiplication, respectively. 
c     Presuming coefficient a(1) is zero (the NPB assumes this,
c     but it is thus not a general case), 3A + 1M may be eliminated,
c     resulting in 12A + 3M.
c     Note that this vectorizes, and is also fine for cache 
c     based machines.  
c-------------------------------------------------------------------*/

    int i3, i2, i1;
    double u1[M], u2[M];
#pragma omp parallel for default(shared) private(i1,i2,i3,u1,u2)
    for (i3 = 1; i3 < n3-1; i3++) {
	for (i2 = 1; i2 < n2-1; i2++) {
            for (i1 = 0; i1 < n1; i1++) {
		u1[i1] = u[i3][i2-1][i1] + u[i3][i2+1][i1]
		       + u[i3-1][i2][i1] + u[i3+1][i2][i1];
		u2[i1] = u[i3-1][i2-1][i1] + u[i3-1][i2+1][i1]
		       + u[i3+1][i2-1][i1] + u[i3+1][i2+1][i1];
	    }
	    for (i1 = 1; i1 < n1-1; i1++) {
		r[i3][i2][i1] = v[i3][i2][i1]
		    - a[0] * u[i3][i2][i1]
/*--------------------------------------------------------------------
c  Assume a(1) = 0      (Enable 2 lines below if a(1) not= 0)
c---------------------------------------------------------------------
c    >                     - a(1) * ( u(i1-1,i2,i3) + u(i1+1,i2,i3)
c    >                              + u1(i1) )
c-------------------------------------------------------------------*/
		- a[2] * ( u2[i1] + u1[i1-1] + u1[i1+1] )
		      - a[3] * ( u2[i1-1] + u2[i1+1] );
	    }
	}
    }

/*--------------------------------------------------------------------
c     exchange boundary data
c--------------------------------------------------------------------*/
    comm3(r,n1,n2,n3,k);

    if (debug_vec[0] >= 1 ) {
	rep_nrm(r,n1,n2,n3,"   resid",k);
    }

    if ( debug_vec[2] >= k ) {
	showall(r,n1,n2,n3);
    }
}