static void zran3(double ***z, int n1, int n2, int n3, int nx, int ny, int k) {
#define MM	10
#define	A	pow(5.0,13)
#define	X	314159265.e0    
    int i0, m0, m1;
    int i1, i2, i3, d1, e1, e2, e3;
    double xx, x0, x1, a1, a2, ai;
    double ten[MM][2], best;
    int i, j1[MM][2], j2[MM][2], j3[MM][2];
    int jg[4][MM][2];
    double rdummy;
    a1 = power( A, nx );
    a2 = power( A, nx*ny );
    zero3(z,n1,n2,n3);
    i = is1-1+nx*(is2-1+ny*(is3-1));
    ai = power( A, i );
    d1 = ie1 - is1 + 1;
    e1 = ie1 - is1 + 2;
    e2 = ie2 - is2 + 2;
    e3 = ie3 - is3 + 2;
    x0 = X;
    rdummy = randlc( &x0, ai );
    for (i3 = 1; i3 < e3; i3++) {
	x1 = x0;
	for (i2 = 1; i2 < e2; i2++) {
            xx = x1;
            vranlc( d1, &xx, A, &(z[i3][i2][0]));
            rdummy = randlc( &x1, a1 );
	}
	rdummy = randlc( &x0, a2 );
    }
    for (i = 0; i < MM; i++) {
	ten[i][1] = 0.0;
	j1[i][1] = 0;
	j2[i][1] = 0;
	j3[i][1] = 0;
	ten[i][0] = 1.0;
	j1[i][0] = 0;
	j2[i][0] = 0;
	j3[i][0] = 0;
    }
    for (i3 = 1; i3 < n3-1; i3++) {
	for (i2 = 1; i2 < n2-1; i2++) {
            for (i1 = 1; i1 < n1-1; i1++) {
		if ( z[i3][i2][i1] > ten[0][1] ) {
		    ten[0][1] = z[i3][i2][i1];
		    j1[0][1] = i1;
		    j2[0][1] = i2;
		    j3[0][1] = i3;
		    bubble( ten, j1, j2, j3, MM, 1 );
		}
		if ( z[i3][i2][i1] < ten[0][0] ) {
		    ten[0][0] = z[i3][i2][i1];
		    j1[0][0] = i1;
		    j2[0][0] = i2;
		    j3[0][0] = i3;
		    bubble( ten, j1, j2, j3, MM, 0 );
		}
	    }
	}
    }
    i1 = MM - 1;
    i0 = MM - 1;
    for (i = MM - 1 ; i >= 0; i--) {
	best = z[j3[i1][1]][j2[i1][1]][j1[i1][1]];
	if (best == z[j3[i1][1]][j2[i1][1]][j1[i1][1]]) {
            jg[0][i][1] = 0;
            jg[1][i][1] = is1 - 1 + j1[i1][1];
            jg[2][i][1] = is2 - 1 + j2[i1][1];
            jg[3][i][1] = is3 - 1 + j3[i1][1];
            i1 = i1-1;
	} else {
            jg[0][i][1] = 0;
            jg[1][i][1] = 0;
            jg[2][i][1] = 0;
            jg[3][i][1] = 0;
	}
	ten[i][1] = best;
	best = z[j3[i0][0]][j2[i0][0]][j1[i0][0]];
	if (best == z[j3[i0][0]][j2[i0][0]][j1[i0][0]]) {
            jg[0][i][0] = 0;
            jg[1][i][0] = is1 - 1 + j1[i0][0];
            jg[2][i][0] = is2 - 1 + j2[i0][0];
            jg[3][i][0] = is3 - 1 + j3[i0][0];
            i0 = i0-1;
	} else {
            jg[0][i][0] = 0;
            jg[1][i][0] = 0;
            jg[2][i][0] = 0;
            jg[3][i][0] = 0;
	}
	ten[i][0] = best;
    }
    m1 = i1+1;
    m0 = i0+1;
    for (i3 = 0; i3 < n3; i3++) {
#pragma omp parallel for private(i2, i1)
	for (i2 = 0; i2 < n2; i2++) {
            for (i1 = 0; i1 < n1; i1++) {
		z[i3][i2][i1] = 0.0;
	    }
	}
    }
    for (i = MM-1; i >= m0; i--) {
	z[j3[i][0]][j2[i][0]][j1[i][0]] = -1.0;
    }
    for (i = MM-1; i >= m1; i--) {
	z[j3[i][1]][j2[i][1]][j1[i][1]] = 1.0;
    } 
    comm3(z,n1,n2,n3,k);
}