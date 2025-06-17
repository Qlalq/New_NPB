#include <omp.h>

static void setbv(void) {
  int i, j, k, iglob, jglob;

  #pragma omp parallel for private(iglob, jglob, i, j) shared(u) 
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (j = 0; j < ny; j++) {
      jglob = j;
      exact(iglob, jglob, 0, &u[i][j][0][0]);
      exact(iglob, jglob, nz-1, &u[i][j][nz-1][0]);
    }
  }

  #pragma omp parallel for private(iglob, k, i) shared(u) 
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (k = 0; k < nz; k++) {
      exact(iglob, 0, k, &u[i][0][k][0]);
    }
  }

  #pragma omp parallel for private(iglob, k, i) shared(u)
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (k = 0; k < nz; k++) {
      exact(iglob, ny0-1, k, &u[i][ny-1][k][0]);
    }
  }

  #pragma omp parallel for private(jglob, k, j) shared(u)
  for (j = 0; j < ny; j++) {
    jglob = j;
    for (k = 0; k < nz; k++) {
      exact(0, jglob, k, &u[0][j][k][0]);
    }
  }

  #pragma omp parallel for private(jglob, k, j) shared(u)
  for (j = 0; j < ny; j++) {
    jglob = j;
    for (k = 0; k < nz; k++) {
      exact(nx0-1, jglob, k, &u[nx-1][j][k][0]);
    }
  }
}