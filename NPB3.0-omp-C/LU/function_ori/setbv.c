static void setbv(void) {

#pragma omp parallel
{

/*--------------------------------------------------------------------
c   set the boundary values of dependent variables
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c   local variables
--------------------------------------------------------------------*/
  int i, j, k;
  int iglob, jglob;

/*--------------------------------------------------------------------
c   set the dependent variable values along the top and bottom faces
--------------------------------------------------------------------*/
#pragma omp for
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (j = 0; j < ny; j++) {
      jglob = j;
      exact( iglob, jglob, 0, &u[i][j][0][0] );
      exact( iglob, jglob, nz-1, &u[i][j][nz-1][0] );
    }
  }

/*--------------------------------------------------------------------
c   set the dependent variable values along north and south faces
--------------------------------------------------------------------*/
#pragma omp for
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (k = 0; k < nz; k++) {
      exact( iglob, 0, k, &u[i][0][k][0] );
    }
  }

#pragma omp for
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (k = 0; k < nz; k++) {
      exact( iglob, ny0-1,  k, &u[i][ny-1][k][0] );
    }
  }

/*--------------------------------------------------------------------
c   set the dependent variable values along east and west faces
--------------------------------------------------------------------*/
#pragma omp for
  for (j = 0; j < ny; j++) {
    jglob = j;
    for (k = 0; k < nz; k++) {
      exact( 0, jglob, k, &u[0][j][k][0] );
    }
  }

#pragma omp for
  for (j = 0; j < ny; j++) {
    jglob = j;
    for (k = 0; k < nz; k++) {
      exact( nx0-1, jglob, k, &u[nx-1][j][k][0] );
    }
  }
}

}