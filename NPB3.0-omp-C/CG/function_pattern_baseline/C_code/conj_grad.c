static void conj_grad (
    int colidx[],    
    int rowstr[],    
    double x[],        
    double z[],        
    double a[],        
    double p[],        
    double q[],        
    double r[],        
    //double w[],        
    double *rnorm )
{
    static int callcount = 0;
    double d, sum, rho, rho0, alpha, beta;
    int i, j, k;
    int cgit, cgitmax = 25;
    rho = 0.0;

    // Initialization of q, z, r, and p vectors
    #pragma omp parallel for
    for (j = 1; j <= naa+1; j++) {
        q[j] = 0.0;
        z[j] = 0.0;
        r[j] = x[j];
        p[j] = r[j];
        //w[j] = 0.0;
    }

    #pragma omp parallel for reduction(+:rho)
    for (j = 1; j <= lastcol-firstcol+1; j++) {
        rho = rho + r[j]*r[j];
    }

    for (cgit = 1; cgit <= cgitmax; cgit++) {
        rho0 = rho;
        d = 0.0;
        rho = 0.0;

        // Calculation of q and the sum for rho
        #pragma omp parallel for private(sum) 
        for (j = 1; j <= lastrow-firstrow+1; j++) {
            sum = 0.0;
            for (k = rowstr[j]; k < rowstr[j+1]; k++) {
                sum = sum + a[k]*p[colidx[k]];
            }
            //w[j] = sum;
            q[j] = sum;
        }

        #pragma omp parallel for reduction(+:d)
        for (j = 1; j <= lastcol-firstcol+1; j++) {
            d = d + p[j]*q[j];
        }

        alpha = rho0 / d;

        #pragma omp parallel for
        for (j = 1; j <= lastcol-firstcol+1; j++) {
            z[j] = z[j] + alpha*p[j];
            r[j] = r[j] - alpha*q[j];
            rho = rho + r[j]*r[j];
        }

        beta = rho / rho0;

        #pragma omp parallel for
        for (j = 1; j <= lastcol-firstcol+1; j++) {
            p[j] = r[j] + beta*p[j];
        }
        callcount++;
    }

    sum = 0.0;

    // Final computation of residual norm
    #pragma omp parallel for reduction(+:sum) private(d)
    for (j = 1; j <= lastrow-firstrow+1; j++) {
        d = 0.0;
        for (k = rowstr[j]; k <= rowstr[j+1]-1; k++) {
            d = d + a[k]*z[colidx[k]];
        }
        r[j] = d;
    }

    #pragma omp parallel for reduction(+:sum) private(d)
    for (j = 1; j <= lastcol-firstcol+1; j++) {
        d = x[j] - r[j];
        sum = sum + d*d;
    }

    (*rnorm) = sqrt(sum);
}