static void norm2u3(double ***r, int n1, int n2, int n3,
		    double *rnm2, double *rnmu, int nx, int ny, int nz) {

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c     norm2u3 evaluates approximations to the L2 norm and the
c     uniform (or L-infinity or Chebyshev) norm, under the
c     assumption that the boundaries are periodic or zero.  Add the
c     boundaries in with half weight (quarter weight on the edges
c     and eighth weight at the corners) for inhomogeneous boundaries.
c-------------------------------------------------------------------*/

    double s = 0.0;
    int i3, i2, i1, n;
    double a = 0.0, tmp = 0.0;

    n = nx*ny*nz;

#pragma omp parallel for default(shared) private(i1,i2,i3,a) reduction(+:s) reduction(max:tmp)
    for (i3 = 1; i3 < n3-1; i3++) {
	for (i2 = 1; i2 < n2-1; i2++) {
            for (i1 = 1; i1 < n1-1; i1++) {
		s = s + r[i3][i2][i1] * r[i3][i2][i1];
		a = fabs(r[i3][i2][i1]);
		if (a > tmp) tmp = a;
	    }
	}
    }
    *rnmu = tmp;
	*rnm2 = sqrt(s/(double)n);
}