static void sprnvc(
    int n,
    int nz,
    double v[],		
    int iv[],		
    int nzloc[],	
    int mark[] ) 	
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
    while (nzv < nz) {
	vecelt = randlc(&tran, amult);
	vecloc = randlc(&tran, amult);
	i = icnvrt(vecloc, nn1) + 1;
	if (i > n) continue;
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