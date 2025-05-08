static void zero3(double ***z, int n1, int n2, int n3) {
    int i1, i2, i3;
#pragma omp parallel for
    for (i3 = 0;i3 < n3; i3++) {
	for (i2 = 0; i2 < n2; i2++) {
            for (i1 = 0; i1 < n1; i1++) {
		z[i3][i2][i1] = 0.0;
	    }
	}
   }
}