static void rprj3( double ***r, int m1k, int m2k, int m3k,
		   double ***s, int m1j, int m2j, int m3j, int k ) {
    int j3, j2, j1, i3, i2, i1, d1, d2, d3;
    double x1[M], y1[M], x2, y2;
    if (m1k == 3) {
        d1 = 2;
    } else {
        d1 = 1;
    }
    if (m2k == 3) {
        d2 = 2;
    } else {
        d2 = 1;
    }
    if (m3k == 3) {
        d3 = 2;
    } else {
        d3 = 1;
    }
    for (j3 = 1; j3 < m3j-1; j3++) {
#pragma omp parallel for collapse(2) private(j1, i3, i2, i1, x2, y2, x1, y1)
	i3 = 2*j3-d3;
	for (j2 = 1; j2 < m2j-1; j2++) {
            i2 = 2*j2-d2;
            for (j1 = 1; j1 < m1j; j1++) {
		i1 = 2*j1-d1;
		x1[i1] = r[i3+1][i2][i1] + r[i3+1][i2+2][i1]
		    + r[i3][i2+1][i1] + r[i3+2][i2+1][i1];
		y1[i1] = r[i3][i2][i1] + r[i3+2][i2][i1]
		    + r[i3][i2+2][i1] + r[i3+2][i2+2][i1];
	    }
            for (j1 = 1; j1 < m1j-1; j1++) {
		i1 = 2*j1-d1;
		y2 = r[i3][i2][i1+1] + r[i3+2][i2][i1+1]
		    + r[i3][i2+2][i1+1] + r[i3+2][i2+2][i1+1];
		x2 = r[i3+1][i2][i1+1] + r[i3+1][i2+2][i1+1]
		    + r[i3][i2+1][i1+1] + r[i3+2][i2+1][i1+1];
		s[j3][j2][j1] =
		    0.5 * r[i3+1][i2+1][i1+1]
		    + 0.25 * ( r[i3+1][i2+1][i1] + r[i3+1][i2+1][i1+2] + x2)
		    + 0.125 * ( x1[i1] + x1[i1+2] + y2)
		    + 0.0625 * ( y1[i1] + y1[i1+2] );
	    }
	}
    }
    comm3(s,m1j,m2j,m3j,k-1);
    if (debug_vec[0] >= 1 ) {
	rep_nrm(s,m1j,m2j,m3j,"   rprj3",k-1);
    }
    if (debug_vec[4] >= k ) {
	showall(s,m1j,m2j,m3j);
    }
}