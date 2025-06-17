static void sprnvc(
    int n,
    int nz,
    double v[],		/* v[1:*] */
    int iv[],		/* iv[1:*] */
    int nzloc[],	/* nzloc[1:n] */
    int mark[] ) 	/* mark[1:n] */
{
    int nn1;
    int nzrow, nzv, ii, i;
    double vecelt, vecloc;

    nzv = 0;
    nzrow = 0;
    nn1 = 1;
    do {
	nn1 = 2 * nn1;
    } while (nn1 < n);

/*--------------------------------------------------------------------
c    nn1 is the smallest power of two not less than n
c-------------------------------------------------------------------*/

    while (nzv < nz) {
	vecelt = randlc(&tran, amult);

/*--------------------------------------------------------------------
c   generate an integer between 1 and n in a portable manner
c-------------------------------------------------------------------*/
	vecloc = randlc(&tran, amult);
	i = icnvrt(vecloc, nn1) + 1;
	if (i > n) continue;

/*--------------------------------------------------------------------
c  was this integer generated already?
c-------------------------------------------------------------------*/
	if (mark[i] == 0) {
	    mark[i] = 1;
	    nzrow = nzrow + 1;
	    nzloc[nzrow] = i;
	    nzv = nzv + 1;
	    v[nzv] = vecelt;
	    iv[nzv] = i;
	}
    }

    for (ii = 1; ii <= nzrow; ii++) {
	i = nzloc[ii];
	mark[i] = 0;
    }
}