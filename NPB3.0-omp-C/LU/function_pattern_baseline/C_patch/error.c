static void error(void) {
  int i, j, k, m;
  int iglob, jglob;
  double  tmp;
  double  u000ijk[5];
  for (m = 0; m < 5; m++) {
    errnm[m] = 0.0;
  }
#pragma omp parallel for \
    private(i, j, k, m, iglob, jglob, tmp, u000ijk) reduction(+:errnm) collapse(3)
  for (i = ist; i <= iend; i++) {
    iglob = i;
    for (j = jst; j <= jend; j++) {
      jglob = j;
      for (k = 1; k <= nz-2; k++) {
	exact( iglob, jglob, k, u000ijk );
	for (m = 0; m < 5; m++) {
	  tmp = ( u000ijk[m] - u[i][j][k][m] );
	  errnm[m] = errnm[m] + tmp *tmp;
	}
      }
    }
  }
  for (m = 0; m < 5; m++) {
    errnm[m] = sqrt ( errnm[m] / ( (nx0-2)*(ny0-2)*(nz0-2) ) );
  }
}