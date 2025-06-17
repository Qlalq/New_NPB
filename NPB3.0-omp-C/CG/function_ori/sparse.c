static void sparse(
    double a[],		/* a[1:*] */
    int colidx[],	/* colidx[1:*] */
    int rowstr[],	/* rowstr[1:*] */
    int n,
    int arow[],		/* arow[1:*] */
    int acol[],		/* acol[1:*] */
    double aelt[],	/* aelt[1:*] */
    int firstrow,
    int lastrow,
    double x[],		/* x[1:n] */
    boolean mark[],	/* mark[1:n] */
    int nzloc[],	/* nzloc[1:n] */
    int nnza)
/*---------------------------------------------------------------------
c       rows range from firstrow to lastrow
c       the rowstr pointers are defined for nrows = lastrow-firstrow+1 values
c---------------------------------------------------------------------*/
{
    int nrows;
    int i, j, jajp1, nza, k, nzrow;
    double xi;

/*--------------------------------------------------------------------
c    how many rows of result
c-------------------------------------------------------------------*/
    nrows = lastrow - firstrow + 1;

/*--------------------------------------------------------------------
c     ...count the number of triples in each row
c-------------------------------------------------------------------*/
#pragma omp parallel for default(shared) private(j)
    for (j = 1; j <= n; j++) {
	rowstr[j] = 0;
	mark[j] = FALSE;
    }
    rowstr[n+1] = 0;
    
    for (nza = 1; nza <= nnza; nza++) {
	j = (arow[nza] - firstrow + 1) + 1;
	rowstr[j] = rowstr[j] + 1;
    }

    rowstr[1] = 1;
    for (j = 2; j <= nrows+1; j++) {
	rowstr[j] = rowstr[j] + rowstr[j-1];
    }

/*---------------------------------------------------------------------
c     ... rowstr(j) now is the location of the first nonzero
c           of row j of a
c---------------------------------------------------------------------*/
    
/*---------------------------------------------------------------------
c     ... preload data pages
c---------------------------------------------------------------------*/
#pragma omp parallel for default(shared) private(k,j)
      for(j = 0;j <= nrows-1;j++) {
         for(k = rowstr[j];k <= rowstr[j+1]-1;k++)
	       a[k] = 0.0;
      }
/*--------------------------------------------------------------------
c     ... do a bucket sort of the triples on the row index
c-------------------------------------------------------------------*/
    for (nza = 1; nza <= nnza; nza++) {
	j = arow[nza] - firstrow + 1;
	k = rowstr[j];
	a[k] = aelt[nza];
	colidx[k] = acol[nza];
	rowstr[j] = rowstr[j] + 1;
    }

/*--------------------------------------------------------------------
c       ... rowstr(j) now points to the first element of row j+1
c-------------------------------------------------------------------*/
    for (j = nrows; j >= 1; j--) {
	rowstr[j+1] = rowstr[j];
    }
    rowstr[1] = 1;

/*--------------------------------------------------------------------
c       ... generate the actual output rows by adding elements
c-------------------------------------------------------------------*/
    nza = 0;
#pragma omp parallel for default(shared) private(i)    
    for (i = 1; i <= n; i++) {
	x[i] = 0.0;
	mark[i] = FALSE;
    }

    jajp1 = rowstr[1];
    for (j = 1; j <= nrows; j++) {
	nzrow = 0;
	
/*--------------------------------------------------------------------
c          ...loop over the jth row of a
c-------------------------------------------------------------------*/
	for (k = jajp1; k < rowstr[j+1]; k++) {
            i = colidx[k];
            x[i] = x[i] + a[k];
            if ( mark[i] == FALSE && x[i] != 0.0) {
		mark[i] = TRUE;
		nzrow = nzrow + 1;
		nzloc[nzrow] = i;
	    }
	}

/*--------------------------------------------------------------------
c          ... extract the nonzeros of this row
c-------------------------------------------------------------------*/
	for (k = 1; k <= nzrow; k++) {
            i = nzloc[k];
            mark[i] = FALSE;
            xi = x[i];
            x[i] = 0.0;
            if (xi != 0.0) {
		nza = nza + 1;
		a[nza] = xi;
		colidx[nza] = i;
	    }
	}
	jajp1 = rowstr[j+1];
	rowstr[j+1] = nza + rowstr[1];
    }
}