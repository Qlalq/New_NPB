static void error(void) {

/*--------------------------------------------------------------------
c
c   compute the solution error
c
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j, k, m;
  int iglob, jglob;
  double  tmp;
  double  u000ijk[5];

  for (m = 0; m < 5; m++) {
    errnm[m] = 0.0;
  }
  
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