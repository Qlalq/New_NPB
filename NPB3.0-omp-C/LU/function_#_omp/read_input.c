static void read_input(void) {
  FILE *fp;
  printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version"
	 " - LU Benchmark\n\n");
  fp = fopen("inputlu.data", "r");
  if (fp != NULL) {
    printf(" Reading from input file inputlu.data\n");
    while(fgetc(fp) != '\n'); while(fgetc(fp) != '\n');
    fscanf(fp, "%d%d", &ipr, &inorm);
    while(fgetc(fp) != '\n');
    while(fgetc(fp) != '\n'); while(fgetc(fp) != '\n');
    fscanf(fp, "%d", &itmax);
    while(fgetc(fp) != '\n');
    while(fgetc(fp) != '\n'); while(fgetc(fp) != '\n');
    fscanf(fp, "%lf", &dt);
    while(fgetc(fp) != '\n');
    while(fgetc(fp) != '\n'); while(fgetc(fp) != '\n');
    fscanf(fp, "%lf", &omega);
    while(fgetc(fp) != '\n');
    while(fgetc(fp) != '\n'); while(fgetc(fp) != '\n');
    fscanf(fp, "%lf%lf%lf%lf%lf",
	   &tolrsd[0], &tolrsd[1], &tolrsd[2], &tolrsd[3], &tolrsd[4]);
    while(fgetc(fp) != '\n');
    while(fgetc(fp) != '\n'); while(fgetc(fp) != '\n');
    fscanf(fp, "%d%d%d", &nx0, &ny0, &nz0);
    while(fgetc(fp) != '\n');
    fclose(fp);
  } else {
    ipr = IPR_DEFAULT;
    inorm = INORM_DEFAULT;
    itmax = ITMAX_DEFAULT;
    dt = DT_DEFAULT;
    omega = OMEGA_DEFAULT;
    tolrsd[0] = TOLRSD1_DEF;
    tolrsd[1] = TOLRSD2_DEF;
    tolrsd[2] = TOLRSD3_DEF;
    tolrsd[3] = TOLRSD4_DEF;
    tolrsd[4] = TOLRSD5_DEF;
    nx0 = ISIZ1;
    ny0 = ISIZ2;
    nz0 = ISIZ3;
  }
  if ( nx0 < 4 || ny0 < 4 || nz0 < 4 ) {
    printf("     PROBLEM SIZE IS TOO SMALL - \n"
	   "     SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5\n");
    exit(1);
  }
  if ( nx0 > ISIZ1 || ny0 > ISIZ2 || nz0 > ISIZ3 ) {
    printf("     PROBLEM SIZE IS TOO LARGE - \n"
	   "     NX, NY AND NZ SHOULD BE EQUAL TO \n"
	   "     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY\n");
    exit(1);
  }
  printf(" Size: %3dx%3dx%3d\n", nx0, ny0, nz0);
  printf(" Iterations: %3d\n", itmax);
}