static void makea(
    int n,
    int nz,
    double a[],		/* a[1:nz] */
    int colidx[],	/* colidx[1:nz] */
    int rowstr[],	/* rowstr[1:n+1] */
    int nonzer,
    int firstrow,
    int lastrow,
    int firstcol,
    int lastcol,
    double rcond,
    int arow[],		/* arow[1:nz] */
    int acol[],		/* acol[1:nz] */
    double aelt[],	/* aelt[1:nz] */
    double v[],		/* v[1:n+1] */
    int iv[],		/* iv[1:2*n+1] */
    double shift )
{
    int i, nnza, iouter, ivelt, ivelt1, irow, nzv;

/*--------------------------------------------------------------------
c      nonzer is approximately  (int(sqrt(nnza /n)));
c-------------------------------------------------------------------*/

    double size, ratio, scale;
    int jcol;

    size = 1.0;
    ratio = pow(rcond, (1.0 / (double)n));
    nnza = 0;

/*---------------------------------------------------------------------
c  Initialize colidx(n+1 .. 2n) to zero.
c  Used by sprnvc to mark nonzero positions
c---------------------------------------------------------------------*/
#pragma omp parallel for default(shared) private(i)
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

/*---------------------------------------------------------------------
c       ... add the identity * rcond to the generated matrix to bound
c           the smallest eigenvalue from below by rcond
c---------------------------------------------------------------------*/
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

/*---------------------------------------------------------------------
c       ... make the sparse matrix from list of elements with duplicates
c           (v and iv are used as  workspace)
c---------------------------------------------------------------------*/
    sparse(a, colidx, rowstr, n, arow, acol, aelt,
	   firstrow, lastrow, v, &(iv[0]), &(iv[n]), nnza);
}