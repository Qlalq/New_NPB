static void norm2u3(double ***r, int n1, int n2, int n3,
		    double *rnm2, double *rnmu, int nx, int ny, int nz) {
    double s = 0.0;
    int i3, i2, i1, n;
    double a = 0.0, tmp = 0.0;
    n = nx*ny*nz;
#pragma omp parallel for private(a) reduction(+:s, max:tmp)
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