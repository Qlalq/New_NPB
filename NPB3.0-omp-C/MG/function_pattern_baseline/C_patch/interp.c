static void interp( double ***z, int mm1, int mm2, int mm3,
		    double ***u, int n1, int n2, int n3, int k ) {
    int i3, i2, i1, d1, d2, d3, t1, t2, t3;
    double z1[M], z2[M], z3[M];
    if ( n1 != 3 && n2 != 3 && n3 != 3 ) {
#pragma omp parallel for default(shared) private(i3, i2, i1, z1, z2, z3) schedule(static)
	for (i3 = 0; i3 < mm3-1; i3++) {
            for (i2 = 0; i2 < mm2-1; i2++) {
		for (i1 = 0; i1 < mm1; i1++) {
		    z1[i1] = z[i3][i2+1][i1] + z[i3][i2][i1];
		    z2[i1] = z[i3+1][i2][i1] + z[i3][i2][i1];
		    z3[i1] = z[i3+1][i2+1][i1] + z[i3+1][i2][i1] + z1[i1];
		}
		for (i1 = 0; i1 < mm1-1; i1++) {
		    u[2*i3][2*i2][2*i1] = u[2*i3][2*i2][2*i1]
			+z[i3][i2][i1];
		    u[2*i3][2*i2][2*i1+1] = u[2*i3][2*i2][2*i1+1]
			+0.5*(z[i3][i2][i1+1]+z[i3][i2][i1]);
		}
		for (i1 = 0; i1 < mm1-1; i1++) {
		    u[2*i3][2*i2+1][2*i1] = u[2*i3][2*i2+1][2*i1]
			+0.5 * z1[i1];
		    u[2*i3][2*i2+1][2*i1+1] = u[2*i3][2*i2+1][2*i1+1]
			+0.25*( z1[i1] + z1[i1+1] );
		}
		for (i1 = 0; i1 < mm1-1; i1++) {
		    u[2*i3+1][2*i2][2*i1] = u[2*i3+1][2*i2][2*i1]
			+0.5 * z2[i1];
		    u[2*i3+1][2*i2][2*i1+1] = u[2*i3+1][2*i2][2*i1+1]
			+0.25*( z2[i1] + z2[i1+1] );
		}
		for (i1 = 0; i1 < mm1-1; i1++) {
		    u[2*i3+1][2*i2+1][2*i1] = u[2*i3+1][2*i2+1][2*i1]
			+0.25* z3[i1];
		    u[2*i3+1][2*i2+1][2*i1+1] = u[2*i3+1][2*i2+1][2*i1+1]
			+0.125*( z3[i1] + z3[i1+1] );
		}
	    }
	}
    } else {
	if (n1 == 3) {
            d1 = 2;
            t1 = 1;
	} else {
            d1 = 1;
            t1 = 0;
	}
	if (n2 == 3) {
            d2 = 2;
            t2 = 1;
	} else {
            d2 = 1;
            t2 = 0;
	}
	if (n3 == 3) {
            d3 = 2;
            t3 = 1;
	} else {
            d3 = 1;
            t3 = 0;
	}
    {
#pragma omp parallel default(shared) private(i3, i2, i1)
#pragma omp for schedule(static) nowait
	for ( i3 = d3; i3 <= mm3-1; i3++) {
            for ( i2 = d2; i2 <= mm2-1; i2++) {
		for ( i1 = d1; i1 <= mm1-1; i1++) { // L64
		    u[2*i3-d3-1][2*i2-d2-1][2*i1-d1-1] =
			u[2*i3-d3-1][2*i2-d2-1][2*i1-d1-1]
			+z[i3-1][i2-1][i1-1];
		}
		for ( i1 = 1; i1 <= mm1-1; i1++) {
		    u[2*i3-d3-1][2*i2-d2-1][2*i1-t1-1] =
			u[2*i3-d3-1][2*i2-d2-1][2*i1-t1-1]
			+0.5*(z[i3-1][i2-1][i1]+z[i3-1][i2-1][i1-1]);
		}
	    }
            for ( i2 = 1; i2 <= mm2-1; i2++) {
		for ( i1 = d1; i1 <= mm1-1; i1++) {
		    u[2*i3-d3-1][2*i2-t2-1][2*i1-d1-1] = // L77
			u[2*i3-d3-1][2*i2-t2-1][2*i1-d1-1]
			+0.5*(z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
		}
		for ( i1 = 1; i1 <= mm1-1; i1++) {
		    u[2*i3-d3-1][2*i2-t2-1][2*i1-t1-1] =
			u[2*i3-d3-1][2*i2-t2-1][2*i1-t1-1]
			+0.25*(z[i3-1][i2][i1]+z[i3-1][i2-1][i1]
			       +z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
		}
	    }
	}
#pragma omp for schedule(static)
	for ( i3 = 1; i3 <= mm3-1; i3++) {

            for ( i2 = d2; i2 <= mm2-1; i2++) {
		for ( i1 = d1; i1 <= mm1-1; i1++) {
		    u[2*i3-t3-1][2*i2-d2-1][2*i1-d1-1] =
			u[2*i3-t3-1][2*i2-d2-1][2*i1-d1-1]
			+0.5*(z[i3][i2-1][i1-1]+z[i3-1][i2-1][i1-1]);
		}
		for ( i1 = 1; i1 <= mm1-1; i1++) {
		    u[2*i3-t3-1][2*i2-d2-1][2*i1-t1-1] =
			u[2*i3-t3-1][2*i2-d2-1][2*i1-t1-1]
			+0.25*(z[i3][i2-1][i1]+z[i3][i2-1][i1-1]
			       +z[i3-1][i2-1][i1]+z[i3-1][i2-1][i1-1]);
		}
	    }
	    for ( i2 = 1; i2 <= mm2-1; i2++) {
		for ( i1 = d1; i1 <= mm1-1; i1++) {
		    u[2*i3-t3-1][2*i2-t2-1][2*i1-d1-1] =
			u[2*i3-t3-1][2*i2-t2-1][2*i1-d1-1]
			+0.25*(z[i3][i2][i1-1]+z[i3][i2-1][i1-1]
			       +z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
		}
		for ( i1 = 1; i1 <= mm1-1; i1++) {
		    u[2*i3-t3-1][2*i2-t2-1][2*i1-t1-1] =
			u[2*i3-t3-1][2*i2-t2-1][2*i1-t1-1]
			+0.125*(z[i3][i2][i1]+z[i3][i2-1][i1]
				+z[i3][i2][i1-1]+z[i3][i2-1][i1-1]
				+z[i3-1][i2][i1]+z[i3-1][i2-1][i1]
				+z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
		}
	    }
	}
    }
    }//end #pragma omp parallel
    if (debug_vec[0] >= 1 ) { // L122
	rep_nrm(z,mm1,mm2,mm3,"z: inter",k-1);
	rep_nrm(u,n1,n2,n3,"u: inter",k);
    }
    if ( debug_vec[5] >= k ) {
	showall(z,mm1,mm2,mm3);
	showall(u,n1,n2,n3);
    }
}