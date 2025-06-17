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

    #pragma omp parallel shared(nzv, nzrow, mark, nzloc, v, iv) private(i, vecelt, vecloc)
    {
        #pragma omp for
        for (i = 0; i < nz; i++) {
            vecelt = randlc(&tran, amult);
            vecloc = randlc(&tran, amult);
            i = icnvrt(vecloc, nn1) + 1;
            if (i > n) continue;
            #pragma omp critical
            {
                if (mark[i] == 0) {
                    mark[i] = 1;
                    nzrow = nzrow + 1;
                    nzloc[nzrow] = i;
                    nzv = nzv + 1;
                    v[nzv] = vecelt;
                    iv[nzv] = i;
                }
            }
        }
        
        #pragma omp for
        for (ii = 1; ii <= nzrow; ii++) {
            i = nzloc[ii];
            mark[i] = 0;
        }
    }
}