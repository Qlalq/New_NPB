static void nonzero(double ***z, int n1, int n2, int n3);
int main(int argc, char *argv[]) {
    int k, it;
    double t, tinit, mflops;
    int nthreads = 1;
    double ****u, ***v, ****r;
    double a[4], c[4];
    double rnm2, rnmu;
    double epsilon = 1.0e-8;
    int n1, n2, n3, nit;
    double verify_value;
    boolean verified;
    int i, j, l;
    FILE *fp;
    timer_clear(T_BENCH);
    timer_clear(T_INIT);
    timer_start(T_INIT);
    printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version"
	   " - MG Benchmark\n\n");
    fp = fopen("mg.input", "r");
    if (fp != NULL) {
	printf(" Reading from input file mg.input\n");
	fscanf(fp, "%d", &lt);
	while(fgetc(fp) != '\n');
	fscanf(fp, "%d%d%d", &nx[lt], &ny[lt], &nz[lt]);
	while(fgetc(fp) != '\n');
	fscanf(fp, "%d", &nit);
	while(fgetc(fp) != '\n');
	for (i = 0; i <= 7; i++) {
	    fscanf(fp, "%d", &debug_vec[i]);
	}
	fclose(fp);
    } else {
	printf(" No input file. Using compiled defaults\n");
	lt = LT_DEFAULT;
	nit = NIT_DEFAULT;
	nx[lt] = NX_DEFAULT;
	ny[lt] = NY_DEFAULT;
	nz[lt] = NZ_DEFAULT;
	for (i = 0; i <= 7; i++) {
	    debug_vec[i] = DEBUG_DEFAULT;
	}
    }
    if ( (nx[lt] != ny[lt]) || (nx[lt] != nz[lt]) ) {
	Class = 'U';
    } else if( nx[lt] == 32 && nit == 4 ) {
	Class = 'S';
    } else if( nx[lt] == 64 && nit == 40 ) {
	Class = 'W';
    } else if( nx[lt] == 256 && nit == 20 ) {
	Class = 'B';
    } else if( nx[lt] == 512 && nit == 20 ) {
	Class = 'C';
    } else if( nx[lt] == 256 && nit == 4 ) {
	Class = 'A';
    } else {
	Class = 'U';
    }
    a[0] = -8.0/3.0;
    a[1] =  0.0;
    a[2] =  1.0/6.0;
    a[3] =  1.0/12.0;
    if (Class == 'A' || Class == 'S' || Class =='W') {
	c[0] =  -3.0/8.0;
	c[1] =  1.0/32.0;
	c[2] =  -1.0/64.0;
	c[3] =   0.0;
    } else {
	c[0] =  -3.0/17.0;
	c[1] =  1.0/33.0;
	c[2] =  -1.0/61.0;
	c[3] =   0.0;
    }
    lb = 1;
    setup(&n1,&n2,&n3,lt);
    u = (double ****)malloc((lt+1)*sizeof(double ***));
    for (l = lt; l >=1; l--) {
	u[l] = (double ***)malloc(m3[l]*sizeof(double **));
	for (k = 0; k < m3[l]; k++) {
	    u[l][k] = (double **)malloc(m2[l]*sizeof(double *));
	    for (j = 0; j < m2[l]; j++) {
		u[l][k][j] = (double *)malloc(m1[l]*sizeof(double));
	    }
	}
    }
    v = (double ***)malloc(m3[lt]*sizeof(double **));
    for (k = 0; k < m3[lt]; k++) {
	v[k] = (double **)malloc(m2[lt]*sizeof(double *));
	for (j = 0; j < m2[lt]; j++) {
	    v[k][j] = (double *)malloc(m1[lt]*sizeof(double));
	}
    }
    r = (double ****)malloc((lt+1)*sizeof(double ***));
    for (l = lt; l >=1; l--) {
	r[l] = (double ***)malloc(m3[l]*sizeof(double **));
	for (k = 0; k < m3[l]; k++) {
	    r[l][k] = (double **)malloc(m2[l]*sizeof(double *));
	    for (j = 0; j < m2[l]; j++) {
		r[l][k][j] = (double *)malloc(m1[l]*sizeof(double));
	    }
	}
    }
    zero3(u[lt],n1,n2,n3);
    zran3(v,n1,n2,n3,nx[lt],ny[lt],lt);
    norm2u3(v,n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);
    printf(" Size: %3dx%3dx%3d (class %1c)\n",
	   nx[lt], ny[lt], nz[lt], Class);
    printf(" Iterations: %3d\n", nit);
    resid(u[lt],v,r[lt],n1,n2,n3,a,lt);
    norm2u3(r[lt],n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);
    mg3P(u,v,r,a,c,n1,n2,n3,lt);
    resid(u[lt],v,r[lt],n1,n2,n3,a,lt);
    setup(&n1,&n2,&n3,lt);
    zero3(u[lt],n1,n2,n3); 
    zran3(v,n1,n2,n3,nx[lt],ny[lt],lt);
    timer_stop(T_INIT);
    timer_start(T_BENCH);
    resid(u[lt],v,r[lt],n1,n2,n3,a,lt);
    norm2u3(r[lt],n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);
    for ( it = 1; it <= nit; it++) {
	mg3P(u,v,r,a,c,n1,n2,n3,lt);
	resid(u[lt],v,r[lt],n1,n2,n3,a,lt);
    }
    norm2u3(r[lt],n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);
{   
#if defined(_OPENMP)
  nthreads = omp_get_num_threads();
#endif 
} 
    timer_stop(T_BENCH);
    t = timer_read(T_BENCH);
    tinit = timer_read(T_INIT);
    verified = FALSE;
    verify_value = 0.0;
    printf(" Initialization time: %15.3f seconds\n", tinit);
    printf(" Benchmark completed\n");
    if (Class != 'U') {
	if (Class == 'S') {
            verify_value = 0.530770700573e-04;
	} else if (Class == 'W') {
            verify_value = 0.250391406439e-17;  
	} else if (Class == 'A') {
            verify_value = 0.2433365309e-5;
        } else if (Class == 'B') {
            verify_value = 0.180056440132e-5;
        } else if (Class == 'C') {
            verify_value = 0.570674826298e-06;
	}
	if ( fabs( rnm2 - verify_value ) <= epsilon ) {
            verified = TRUE;
	    printf(" VERIFICATION SUCCESSFUL\n");
	    printf(" L2 Norm is %20.12e\n", rnm2);
	    printf(" Error is   %20.12e\n", rnm2 - verify_value);
	} else {
            verified = FALSE;
	    printf(" VERIFICATION FAILED\n");
	    printf(" L2 Norm is             %20.12e\n", rnm2);
	    printf(" The correct L2 Norm is %20.12e\n", verify_value);
	}
    } else {
	verified = FALSE;
	printf(" Problem size unknown\n");
	printf(" NO VERIFICATION PERFORMED\n");
    }
    if ( t != 0.0 ) {
	int nn = nx[lt]*ny[lt]*nz[lt];
	mflops = 58.*nit*nn*1.0e-6 / t;
    } else {
	mflops = 0.0;
    }
    c_print_results("MG", Class, nx[lt], ny[lt], nz[lt], 
		    nit, nthreads, t, mflops, "          floating point", 
		    verified, NPBVERSION, COMPILETIME,
		    CS1, CS2, CS3, CS4, CS5, CS6, CS7);
}