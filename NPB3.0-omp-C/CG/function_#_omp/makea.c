static void makea(
    int n,
    int nz,
    double a[],		
    int colidx[],	
    int rowstr[],	
    int nonzer,
    int firstrow,
    int lastrow,
    int firstcol,
    int lastcol,
    double rcond,
    int arow[],		
    int acol[],		
    double aelt[],	
    double v[],		
    int iv[],		
    double shift )
{
    int i, nnza, iouter, ivelt, ivelt1, irow, nzv;
    double size, ratio, scale;
    int jcol;
    size = 1.0;
    ratio = pow(rcond, (1.0 / (double)n));
    nnza = 0;
    for (i = 1; i <= n; i++) {
	colidx[n+i] = 0;
    }
    for (iouter = 1; iouter <= n; iouter++) {
	nzv = nonzer;
	sprnvc(n, nzv, v, iv, &(colidx[0]), &(colidx[n]));
	vecset(n, v, iv, &nzv, iouter, 0.5);
	for (ivelt = 1; ivelt <= nzv; ivelt++) {
	    jcol = iv[ivelt];
	    if (jcol >= firstcol && jcol <= lastcol) {
		scale = size * v[ivelt];
		for (ivelt1 = 1; ivelt1 <= nzv; ivelt1++) {
	            irow = iv[ivelt1];
                    if (irow >= firstrow && irow <= lastrow) {
			nnza = nnza + 1;
			if (nnza > nz) {
			    printf("Space for matrix elements exceeded in"
				   " makea\n");
			    printf("nnza, nzmax = %d, %d\n", nnza, nz);
			    printf("iouter = %d\n", iouter);
			    exit(1);
			}
			acol[nnza] = jcol;
			arow[nnza] = irow;
			aelt[nnza] = v[ivelt1] * scale;
		    }
		}
	    }
	}
	size = size * ratio;
    }
    for (i = firstrow; i <= lastrow; i++) {
	if (i >= firstcol && i <= lastcol) {
	    iouter = n + i;
	    nnza = nnza + 1;
	    if (nnza > nz) {
		printf("Space for matrix elements exceeded in makea\n");
		printf("nnza, nzmax = %d, %d\n", nnza, nz);
		printf("iouter = %d\n", iouter);
		exit(1);
	    }
	    acol[nnza] = i;
	    arow[nnza] = i;
	    aelt[nnza] = rcond - shift;
	}
    }
    sparse(a, colidx, rowstr, n, arow, acol, aelt,
	   firstrow, lastrow, v, &(iv[0]), &(iv[n]), nnza);
}