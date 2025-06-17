static void psinv( double ***r, double ***u, int n1, int n2, int n3,
		   double c[4], int k) {

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c     psinv applies an approximate inverse as smoother:  u = u + Cr
c
c     This  implementation costs  15A + 4M per result, where
c     A and M denote the costs of Addition and Multiplication.  
c     Presuming coefficient c(3) is zero (the NPB assumes this,
c     but it is thus not a general case), 2A + 1M may be eliminated,
c     resulting in 13A + 3M.
c     Note that this vectorizes, and is also fine for cache 
c     based machines.  
c-------------------------------------------------------------------*/

    int i3, i2, i1;
    double r1[M], r2[M];
#pragma omp parallel for default(shared) private(i1,i2,i3,r1,r2)   
    for (i3 = 1; i3 < n3-1; i3++) {
	for (i2 = 1; i2 < n2-1; i2++) {
            for (i1 = 0; i1 < n1; i1++) {
		r1[i1] = r[i3][i2-1][i1] + r[i3][i2+1][i1]
		    + r[i3-1][i2][i1] + r[i3+1][i2][i1];
		r2[i1] = r[i3-1][i2-1][i1] + r[i3-1][i2+1][i1]
		    + r[i3+1][i2-1][i1] + r[i3+1][i2+1][i1];
	    }
            for (i1 = 1; i1 < n1-1; i1++) {
		u[i3][i2][i1] = u[i3][i2][i1]
		    + c[0] * r[i3][i2][i1]
		    + c[1] * ( r[i3][i2][i1-1] + r[i3][i2][i1+1]
			       + r1[i1] )
		    + c[2] * ( r2[i1] + r1[i1-1] + r1[i1+1] );
/*--------------------------------------------------------------------
c  Assume c(3) = 0    (Enable line below if c(3) not= 0)
c---------------------------------------------------------------------
c    >                     + c(3) * ( r2(i1-1) + r2(i1+1) )
c-------------------------------------------------------------------*/
	    }
	}
    }

/*--------------------------------------------------------------------
c     exchange boundary points
c-------------------------------------------------------------------*/
    comm3(u,n1,n2,n3,k);

    if (debug_vec[0] >= 1 ) {
	rep_nrm(u,n1,n2,n3,"   psinv",k);
    }

    if ( debug_vec[3] >= k ) {
	showall(u,n1,n2,n3);
    }
}