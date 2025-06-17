static void sparse(
    double a[],		
    int colidx[],	
    int rowstr[],	
    int n,
    int arow[],		
    int acol[],		
    double aelt[],	
    int firstrow,
    int lastrow,
    double x[],		
    boolean mark[],	
    int nzloc[],	
    int nnza)
{
    int nrows;
    int i, j, jajp1, nza, k, nzrow;
    double xi;
    nrows = lastrow - firstrow + 1;
    
    #pragma omp parallel for
    for (j = 1; j <= n; j++) {
        rowstr[j] = 0;
        mark[j] = FALSE;
    }

    rowstr[n+1] = 0;
    
    #pragma omp parallel for private(j)
    for (nza = 1; nza <= nnza; nza++) {
        j = (arow[nza] - firstrow + 1) + 1;
        #pragma omp atomic
        rowstr[j] = rowstr[j] + 1;
    }

    rowstr[1] = 1;
    
    #pragma omp parallel for
    for (j = 2; j <= nrows+1; j++) {
        rowstr[j] = rowstr[j] + rowstr[j-1];
    }
    
    #pragma omp parallel for private(k)
    for(j = 0; j <= nrows-1; j++) {
        for(k = rowstr[j]; k <= rowstr[j+1]-1; k++) {
            a[k] = 0.0;
        }
    }
    
    #pragma omp parallel for private(j, k)
    for (nza = 1; nza <= nnza; nza++) {
        j = arow[nza] - firstrow + 1;
        k = rowstr[j];
        a[k] = aelt[nza];
        colidx[k] = acol[nza];
        #pragma omp atomic
        rowstr[j] = rowstr[j] + 1;
    }

    #pragma omp parallel for
    for (j = nrows; j >= 1; j--) {
        rowstr[j+1] = rowstr[j];
    }
    
    rowstr[1] = 1;
    nza = 0;
    
    #pragma omp parallel for
    for (i = 1; i <= n; i++) {
        x[i] = 0.0;
        mark[i] = FALSE;
    }

    jajp1 = rowstr[1];
    
    #pragma omp parallel for private(k, nzrow, i, xi)
    for (j = 1; j <= nrows; j++) {
        nzrow = 0;
        for (k = jajp1; k < rowstr[j+1]; k++) {
            i = colidx[k];
            x[i] = x[i] + a[k];
            if ( mark[i] == FALSE && x[i] != 0.0) {
                mark[i] = TRUE;
                nzrow = nzrow + 1;
                nzloc[nzrow] = i;
            }
        }
        
        for (k = 1; k <= nzrow; k++) {
            i = nzloc[k];
            mark[i] = FALSE;
            xi = x[i];
            x[i] = 0.0;
            if (xi != 0.0) {
                #pragma omp atomic
                nza = nza + 1;
                a[nza] = xi;
                colidx[nza] = i;
            }
        }
        jajp1 = rowstr[j+1];
        rowstr[j+1] = nza + rowstr[1];
    }
}