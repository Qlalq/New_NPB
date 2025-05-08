static void domain(void) {
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