static void psinv( double ***r, double ***u, int n1, int n2, int n3,
		   double c[4], int k) {
    int i3, i2, i1;
    double r1[M], r2[M];
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
	    }
	}
    }
    comm3(u,n1,n2,n3,k);
    if (debug_vec[0] >= 1 ) {
	rep_nrm(u,n1,n2,n3,"   psinv",k);
    }
    if ( debug_vec[3] >= k ) {
	showall(u,n1,n2,n3);
    }
}