static void evolve(dcomplex u0[NZ][NY][NX], dcomplex u1[NZ][NY][NX],
		   int t, int indexmap[NZ][NY][NX], int d[3]) {
    int i, j, k;
    for (k = 0; k < d[2]; k++) {
	for (j = 0; j < d[1]; j++) {
            for (i = 0; i < d[0]; i++) {
	      crmul(u1[k][j][i], u0[k][j][i], ex[t*indexmap[k][j][i]]);
	    }
	}
    }
}