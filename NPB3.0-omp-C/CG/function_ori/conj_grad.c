static void conj_grad (
    int colidx[],	/* colidx[1:nzz] */
    int rowstr[],	/* rowstr[1:naa+1] */
    double x[],		/* x[*] */
    double z[],		/* z[*] */
    double a[],		/* a[1:nzz] */
    double p[],		/* p[*] */
    double q[],		/* q[*] */
    double r[],		/* r[*] */
    //double w[],		/* w[*] */
    double *rnorm )
/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/
    
/*---------------------------------------------------------------------
c  Floaging point arrays here are named as in NPB1 spec discussion of 
c  CG algorithm
c---------------------------------------------------------------------*/
{
    static int callcount = 0;
    double d, sum, rho, rho0, alpha, beta;
    int i, j, k;
    int cgit, cgitmax = 25;

    rho = 0.0;
#pragma omp parallel default(shared) private(j,sum) shared(rho,naa)
    
/*--------------------------------------------------------------------
c  Initialize the CG algorithm:
c-------------------------------------------------------------------*/
{
#pragma omp for
    for (j = 1; j <= naa+1; j++) {
	q[j] = 0.0;
	z[j] = 0.0;
	r[j] = x[j];
	p[j] = r[j];
	//w[j] = 0.0;
    }

/*--------------------------------------------------------------------
c  rho = r.r
c  Now, obtain the norm of r: First, sum squares of r elements locally...
c-------------------------------------------------------------------*/
#pragma omp for reduction(+:rho)
    for (j = 1; j <= lastcol-firstcol+1; j++) {
	rho = rho + r[j]*r[j];
    }
}/* end omp parallel */
/*--------------------------------------------------------------------
c---->
c  The conj grad iteration loop
c---->
c-------------------------------------------------------------------*/
    for (cgit = 1; cgit <= cgitmax; cgit++) {
      rho0 = rho;
      d = 0.0;
      rho = 0.0;
#pragma omp parallel default(shared) private(j,k,sum,alpha,beta) shared(d,rho0,rho)
{
      
/*--------------------------------------------------------------------
c  q = A.p
c  The partition submatrix-vector multiply: use workspace w
c---------------------------------------------------------------------
C
C  NOTE: this version of the multiply is actually (slightly: maybe %5) 
C        faster on the sp2 on 16 nodes than is the unrolled-by-2 version 
C        below.   On the Cray t3d, the reverse is true, i.e., the 
C        unrolled-by-two version is some 10% faster.  
C        The unrolled-by-8 version below is significantly faster
C        on the Cray t3d - overall speed of code is 1.5 times faster.
*/

/* rolled version */    
#pragma omp for 
	for (j = 1; j <= lastrow-firstrow+1; j++) {
            sum = 0.0;
	    for (k = rowstr[j]; k < rowstr[j+1]; k++) {
		sum = sum + a[k]*p[colidx[k]];
	    }
            //w[j] = sum;
            q[j] = sum;
	}
	
/* unrolled-by-two version
#pragma omp for private(i,k)
        for (j = 1; j <= lastrow-firstrow+1; j++) {
	    int iresidue;
	    double sum1, sum2;
	    i = rowstr[j]; 
            iresidue = (rowstr[j+1]-i) % 2;
            sum1 = 0.0;
            sum2 = 0.0;
            if (iresidue == 1) sum1 = sum1 + a[i]*p[colidx[i]];
	    for (k = i+iresidue; k <= rowstr[j+1]-2; k += 2) {
		sum1 = sum1 + a[k]   * p[colidx[k]];
		sum2 = sum2 + a[k+1] * p[colidx[k+1]];
	    }
            w[j] = sum1 + sum2;
        }
*/
/* unrolled-by-8 version
#pragma omp for private(i,k,sum)
        for (j = 1; j <= lastrow-firstrow+1; j++) {
	    int iresidue;
            i = rowstr[j]; 
            iresidue = (rowstr[j+1]-i) % 8;
            sum = 0.0;
            for (k = i; k <= i+iresidue-1; k++) {
                sum = sum +  a[k] * p[colidx[k]];
            }
            for (k = i+iresidue; k <= rowstr[j+1]-8; k += 8) {
                sum = sum + a[k  ] * p[colidx[k  ]]
                          + a[k+1] * p[colidx[k+1]]
                          + a[k+2] * p[colidx[k+2]]
                          + a[k+3] * p[colidx[k+3]]
                          + a[k+4] * p[colidx[k+4]]
                          + a[k+5] * p[colidx[k+5]]
                          + a[k+6] * p[colidx[k+6]]
                          + a[k+7] * p[colidx[k+7]];
            }
            w[j] = sum;
        }
*/
/*	
#pragma omp for
	for (j = 1; j <= lastcol-firstcol+1; j++) {
            q[j] = w[j];
	}
*/
/*--------------------------------------------------------------------
c  Clear w for reuse...
c-------------------------------------------------------------------*/
/*
#pragma omp for	nowait
	for (j = 1; j <= lastcol-firstcol+1; j++) {
            w[j] = 0.0;
	}
*/
/*--------------------------------------------------------------------
c  Obtain p.q
c-------------------------------------------------------------------*/
#pragma omp for reduction(+:d)
	for (j = 1; j <= lastcol-firstcol+1; j++) {
            d = d + p[j]*q[j];
	}
#pragma omp barrier
/*--------------------------------------------------------------------
c  Obtain alpha = rho / (p.q)
c-------------------------------------------------------------------*/
//#pragma omp single	
	alpha = rho0 / d;

/*--------------------------------------------------------------------
c  Save a temporary of rho
c-------------------------------------------------------------------*/
	/*	rho0 = rho;*/

/*---------------------------------------------------------------------
c  Obtain z = z + alpha*p
c  and    r = r - alpha*q
c---------------------------------------------------------------------*/
#pragma omp for reduction(+:rho)	
	for (j = 1; j <= lastcol-firstcol+1; j++) {
            z[j] = z[j] + alpha*p[j];
            r[j] = r[j] - alpha*q[j];
//	}
            
/*---------------------------------------------------------------------
c  rho = r.r
c  Now, obtain the norm of r: First, sum squares of r elements locally...
c---------------------------------------------------------------------*/
/*
#pragma omp for
	for (j = 1; j <= lastcol-firstcol+1; j++) {*/
            rho = rho + r[j]*r[j];
	}
//#pragma omp barrier

/*--------------------------------------------------------------------
c  Obtain beta:
c-------------------------------------------------------------------*/
//#pragma omp single	
	beta = rho / rho0;

/*--------------------------------------------------------------------
c  p = r + beta*p
c-------------------------------------------------------------------*/
#pragma omp for nowait
	for (j = 1; j <= lastcol-firstcol+1; j++) {
            p[j] = r[j] + beta*p[j];
	}
    callcount++;
    } /* end omp parallel */
    } /* end of do cgit=1,cgitmax */

/*---------------------------------------------------------------------
c  Compute residual norm explicitly:  ||r|| = ||x - A.z||
c  First, form A.z
c  The partition submatrix-vector multiply
c---------------------------------------------------------------------*/
    sum = 0.0;
    
#pragma omp parallel default(shared) private(j,d) shared(sum)
{
#pragma omp for //private(d, k)
    for (j = 1; j <= lastrow-firstrow+1; j++) {
	d = 0.0;
	for (k = rowstr[j]; k <= rowstr[j+1]-1; k++) {
            d = d + a[k]*z[colidx[k]];
	}
	r[j] = d;
    }

/*--------------------------------------------------------------------
c  At this point, r contains A.z
c-------------------------------------------------------------------*/
#pragma omp for reduction(+:sum)
    for (j = 1; j <= lastcol-firstcol+1; j++) {
	d = x[j] - r[j];
	sum = sum + d*d;
    }
} //end omp parallel
    (*rnorm) = sqrt(sum);
}