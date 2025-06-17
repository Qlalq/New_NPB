static void setbv(void) {
{
  int i, j, k;
  int iglob, jglob;
#pragma omp parallel for private(i, j, iglob, jglob) schedule(static)
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (j = 0; j < ny; j++) {
      jglob = j;
      exact( iglob, jglob, 0, &u[i][j][0][0] );
      exact( iglob, jglob, nz-1, &u[i][j][nz-1][0] );
    }
  }
#pragma omp parallel for private(i, k, iglob) schedule(static)
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (k = 0; k < nz; k++) {
      exact( iglob, 0, k, &u[i][0][k][0] );
    }
  }
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (k = 0; k < nz; k++) {
      exact( iglob, ny0-1,  k, &u[i][ny-1][k][0] );
    }
  }
  for (j = 0; j < ny; j++) {
    jglob = j;
    for (k = 0; k < nz; k++) {
      exact( 0, jglob, k, &u[0][j][k][0] );
    }
  }
  for (j = 0; j < ny; j++) {
    jglob = j;
    for (k = 0; k < nz; k++) {
      exact( nx0-1, jglob, k, &u[nx-1][j][k][0] );
    }
  }
}
}
static void setbv(void) {
{
  int i, j, k;
  int iglob, jglob;
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (j = 0; j < ny; j++) {
      jglob = j;
      exact( iglob, jglob, 0, &u[i][j][0][0] );
      exact( iglob, jglob, nz-1, &u[i][j][nz-1][0] );
    }
  }
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (k = 0; k < nz; k++) {
      exact( iglob, 0, k, &u[i][0][k][0] );
    }
  }
#pragma omp parallel for private(i, k, iglob) schedule(static)
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (k = 0; k < nz; k++) {
      exact( iglob, ny0-1,  k, &u[i][ny-1][k][0] );
    }
  }
#pragma omp parallel for private(j, k, jglob) schedule(static)
  for (j = 0; j < ny; j++) {
    jglob = j;
    for (k = 0; k < nz; k++) {
      exact( 0, jglob, k, &u[0][j][k][0] );
    }
  }
#pragma omp parallel for private(j, k, jglob) schedule(static)
  for (j = 0; j < ny; j++) {
    jglob = j;
    for (k = 0; k < nz; k++) {
      exact( nx0-1, jglob, k, &u[nx-1][j][k][0] );
    }
  }
}
}