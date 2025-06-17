static void resid( double ***u, double ***v, double ***r,
		   int n1, int n2, int n3, double a[4], int k ) {
    int i3, i2, i1;
    double u1[M], u2[M];
#pragma omp parallel for private(i3, i2, i1, u1, u2)
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
		- a[2] * ( u2[i1] + u1[i1-1] + u1[i1+1] )
		      - a[3] * ( u2[i1-1] + u2[i1+1] );
	    }
	}
    }
    comm3(r,n1,n2,n3,k);
    if (debug_vec[0] >= 1 ) {
	rep_nrm(r,n1,n2,n3,"   resid",k);
    }
    if ( debug_vec[2] >= k ) {
	showall(r,n1,n2,n3);
    }
}