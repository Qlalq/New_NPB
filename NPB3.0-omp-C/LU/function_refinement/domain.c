static void domain(void) {
  // This function performs initialization and boundary checks.
  // It does not contain significant computational loops or data dependencies
  // that would typically be targeted for parallelization using OpenMP primitives
  // like parallel for, reduction, etc.
  // The variables modified (nx, ny, nz, ist, iend, jst, jend) are likely setup
  // parameters used *by* subsequent parallel regions.
  // Therefore, this function is inherently sequential and doesn't require
  // refactoring according to the principles of decoupling data dependencies,
  // load balancing, or reducing synchronization overhead for OpenMP within *this* function.
  // It is suitable for execution by a single thread before entering parallel regions.

  nx = nx0;
  ny = ny0;
  nz = nz0;

  if ( nx < 4 || ny < 4 || nz < 4 ) {
    printf("     SUBDOMAIN SIZE IS TOO SMALL - \n"
	   "     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n"
	   "     SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL\n"
	   "     TO 4 THEY ARE CURRENTLY%3d%3d%3d\n", nx, ny, nz);
    exit(1);
  }

  if ( nx > ISIZ1 || ny > ISIZ2 || nz > ISIZ3 ) {
    printf("     SUBDOMAIN SIZE IS TOO LARGE - \n"
	   "     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n"
	   "     SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL TO \n"
	   "     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY.  THEY ARE\n"
	   "     CURRENTLY%4d%4d%4d\n", nx, ny, nz);
    exit(1);
  }

  ist = 1;
  iend = nx - 2;
  jst = 1;
  jend = ny - 2;
}