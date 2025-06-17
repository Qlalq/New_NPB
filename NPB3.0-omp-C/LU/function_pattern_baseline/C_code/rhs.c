static void rhs(void) {
  int i, j, k, m;
  int L1, L2;
  int ist1, iend1;
  int jst1, jend1;
  double  q;
  double  u21, u31, u41;
  double  tmp;
  double  u21i, u31i, u41i, u51i;
  double  u21j, u31j, u41j, u51j;
  double  u21k, u31k, u41k, u51k;
  double  u21im1, u31im1, u41im1, u51im1;
  double  u21jm1, u31jm1, u41jm1, u51jm1;
  double  u21km1, u31km1, u41km1, u51km1;

  /*-------------------------------------------------------------------------
    zero out residuals: rsd = -frct
  -------------------------------------------------------------------------*/
#pragma omp parallel for collapse(3) private(i,j,k,m) schedule(static)
  for (i = 0; i <= nx-1; i++) {
    for (j = 0; j <= ny-1; j++) {
      for (k = 0; k <= nz-1; k++) {
        for (m = 0; m < 5; m++) {
          rsd[i][j][k][m] = - frct[i][j][k][m];
        }
      }
    }
  }

  /*-------------------------------------------------------------------------
    compute fluxes in X-direction
  -------------------------------------------------------------------------*/
  L1 = 0;
  L2 = nx-1;
#pragma omp parallel for collapse(3) private(i,j,k,u21,q) schedule(static)
  for (i = L1; i <= L2; i++) {
    for (j = jst; j <= jend; j++) {
      for (k = 1; k <= nz - 2; k++) {
        flux[i][j][k][0] = u[i][j][k][1];
        u21 = u[i][j][k][1] / u[i][j][k][0];
        q = 0.5 * (u[i][j][k][1]*u[i][j][k][1]
                + u[i][j][k][2]*u[i][j][k][2]
                + u[i][j][k][3]*u[i][j][k][3])
            / u[i][j][k][0];
        flux[i][j][k][1] = u[i][j][k][1]*u21 + C2*(u[i][j][k][4] - q);
        flux[i][j][k][2] = u[i][j][k][2]*u21;
        flux[i][j][k][3] = u[i][j][k][3]*u21;
        flux[i][j][k][4] = (C1*u[i][j][k][4] - C2*q)*u21;
      }
    }
  }

  /*-------------------------------------------------------------------------
    add X‐flux differences to residuals and apply X‐viscous terms
  -------------------------------------------------------------------------*/
#pragma omp parallel for collapse(3) private(i,j,k,m,tmp, \
      u21i,u31i,u41i,u51i, u21im1,u31im1,u41im1,u51im1) schedule(static)
  for (j = jst; j <= jend; j++) {
    for (k = 1; k <= nz - 2; k++) {
      /* flux difference */
      for (i = ist; i <= iend; i++) {
        for (m = 0; m < 5; m++) {
          rsd[i][j][k][m] -= tx2 * (flux[i+1][j][k][m] - flux[i-1][j][k][m]);
        }
      }
      /* viscous flux */
      L2 = nx-1;
      for (i = ist; i <= L2; i++) {
        tmp = 1.0 / u[i][j][k][0];
        u21i = tmp * u[i][j][k][1];
        u31i = tmp * u[i][j][k][2];
        u41i = tmp * u[i][j][k][3];
        u51i = tmp * u[i][j][k][4];
        tmp = 1.0 / u[i-1][j][k][0];
        u21im1 = tmp * u[i-1][j][k][1];
        u31im1 = tmp * u[i-1][j][k][2];
        u41im1 = tmp * u[i-1][j][k][3];
        u51im1 = tmp * u[i-1][j][k][4];
        flux[i][j][k][1] = (4.0/3.0)*tx3*(u21i-u21im1);
        flux[i][j][k][2] = tx3*(u31i-u31im1);
        flux[i][j][k][3] = tx3*(u41i-u41im1);
        flux[i][j][k][4] = 0.5*(1.0-C1*C5)*tx3*
                           ((u21i*u21i+u31i*u31i+u41i*u41i)
                          -(u21im1*u21im1+u31im1*u31im1+u41im1*u41im1))
                          + (1.0/6.0)*tx3*(u21i*u21i - u21im1*u21im1)
                          + C1*C5*tx3*(u51i-u51im1);
      }
      /* add viscous terms and second‐derivative DS */
      for (i = ist; i <= iend; i++) {
        rsd[i][j][k][0] += dx1*tx1*(u[i-1][j][k][0]
                                 -2.0*u[i][j][k][0]
                                  +u[i+1][j][k][0]);
        rsd[i][j][k][1] += tx3*C3*C4*(flux[i+1][j][k][1]-flux[i][j][k][1])
                        +dx2*tx1*(u[i-1][j][k][1]
                                 -2.0*u[i][j][k][1]
                                  +u[i+1][j][k][1]);
        rsd[i][j][k][2] += tx3*C3*C4*(flux[i+1][j][k][2]-flux[i][j][k][2])
                        +dx3*tx1*(u[i-1][j][k][2]
                                 -2.0*u[i][j][k][2]
                                  +u[i+1][j][k][2]);
        rsd[i][j][k][3] += tx3*C3*C4*(flux[i+1][j][k][3]-flux[i][j][k][3])
                        +dx4*tx1*(u[i-1][j][k][3]
                                 -2.0*u[i][j][k][3]
                                  +u[i+1][j][k][3]);
        rsd[i][j][k][4] += tx3*C3*C4*(flux[i+1][j][k][4]-flux[i][j][k][4])
                        +dx5*tx1*(u[i-1][j][k][4]
                                 -2.0*u[i][j][k][4]
                                  +u[i+1][j][k][4]);
      }
      /* boundary DS (i=1,2) */
      for (m = 0; m < 5; m++) {
        rsd[1][j][k][m]   -= dssp*(+5.0*u[1][j][k][m]
                                 -4.0*u[2][j][k][m]
                                  +   u[3][j][k][m]);
        rsd[2][j][k][m]   -= dssp*(-4.0*u[1][j][k][m]
                                 +6.0*u[2][j][k][m]
                                 -4.0*u[3][j][k][m]
                                  +   u[4][j][k][m]);
      }
      /* internal DS */
      ist1 = 3; iend1 = nx-4;
      for (i = ist1; i <= iend1; i++)
        for (m = 0; m < 5; m++)
          rsd[i][j][k][m] -= dssp*(u[i-2][j][k][m]
                                -4.0*u[i-1][j][k][m]
                                +6.0*u[i][j][k][m]
                                -4.0*u[i+1][j][k][m]
                                 +u[i+2][j][k][m]);
      /* boundary DS (i=nx-3,nx-2) */
      for (m = 0; m < 5; m++) {
        rsd[nx-3][j][k][m] -= dssp*( u[nx-5][j][k][m]
                                  -4.0*u[nx-4][j][k][m]
                                  +6.0*u[nx-3][j][k][m]
                                  -4.0*u[nx-2][j][k][m]);
        rsd[nx-2][j][k][m] -= dssp*( u[nx-4][j][k][m]
                                  -4.0*u[nx-3][j][k][m]
                                  +5.0*u[nx-2][j][k][m]);
      }
    }
  }

  /*-------------------------------------------------------------------------
    compute fluxes in Y‐direction
  -------------------------------------------------------------------------*/
  L1 = 0;
  L2 = ny-1;
#pragma omp parallel for collapse(3) private(i,j,k,u31,q) schedule(static)
  for (i = ist; i <= iend; i++) {
    for (j = L1; j <= L2; j++) {
      for (k = 1; k <= nz - 2; k++) {
        flux[i][j][k][0] = u[i][j][k][2];
        u31 = u[i][j][k][2] / u[i][j][k][0];
        q = 0.5 * (u[i][j][k][1]*u[i][j][k][1]
                + u[i][j][k][2]*u[i][j][k][2]
                + u[i][j][k][3]*u[i][j][k][3])
            / u[i][j][k][0];
        flux[i][j][k][1] = u[i][j][k][1]*u31;
        flux[i][j][k][2] = u[i][j][k][2]*u31 + C2*(u[i][j][k][4]-q);
        flux[i][j][k][3] = u[i][j][k][3]*u31;
        flux[i][j][k][4] = (C1*u[i][j][k][4] - C2*q)*u31;
      }
    }
  }

  /*-------------------------------------------------------------------------
    add Y‐flux differences and Y‐viscous terms
  -------------------------------------------------------------------------*/
#pragma omp parallel for collapse(3) private(i,j,k,m,tmp, \
      u21j,u31j,u41j,u51j, u21jm1,u31jm1,u41jm1,u51jm1) schedule(static)
  for (i = ist; i <= iend; i++) {
    for (k = 1; k <= nz - 2; k++) {
      /* flux diff */
      for (j = jst; j <= jend; j++) {
        for (m = 0; m < 5; m++) {
          rsd[i][j][k][m] -= ty2*(flux[i][j+1][k][m] - flux[i][j-1][k][m]);
        }
      }
      /* viscous flux */
      L2 = ny-1;
      for (j = jst; j <= L2; j++) {
        tmp = 1.0 / u[i][j][k][0];
        u21j = tmp * u[i][j][k][1];
        u31j = tmp * u[i][j][k][2];
        u41j = tmp * u[i][j][k][3];
        u51j = tmp * u[i][j][k][4];
        tmp = 1.0 / u[i][j-1][k][0];
        u21jm1 = tmp * u[i][j-1][k][1];
        u31jm1 = tmp * u[i][j-1][k][2];
        u41jm1 = tmp * u[i][j-1][k][3];
        u51jm1 = tmp * u[i][j-1][k][4];
        flux[i][j][k][1] = ty3*(u21j - u21jm1);
        flux[i][j][k][2] = (4.0/3.0)*ty3*(u31j-u31jm1);
        flux[i][j][k][3] = ty3*(u41j - u41jm1);
        flux[i][j][k][4] = 0.5*(1.0-C1*C5)*ty3*
                           ((u21j*u21j+u31j*u31j+u41j*u41j)
                          -(u21jm1*u21jm1+u31jm1*u31jm1+u41jm1*u41jm1))
                          + (1.0/6.0)*ty3*(u31j*u31j - u31jm1*u31jm1)
                          + C1*C5*ty3*(u51j-u51jm1);
      }
      /* add viscous terms and DS in Y */
      for (j = jst; j <= jend; j++) {
        rsd[i][j][k][0] += dy1*ty1*(u[i][j-1][k][0]
                                 -2.0*u[i][j][k][0]
                                  +u[i][j+1][k][0]);
        rsd[i][j][k][1] += ty3*C3*C4*(flux[i][j+1][k][1]-flux[i][j][k][1])
                        +dy2*ty1*(u[i][j-1][k][1]
                                 -2.0*u[i][j][k][1]
                                  +u[i][j+1][k][1]);
        rsd[i][j][k][2] += ty3*C3*C4*(flux[i][j+1][k][2]-flux[i][j][k][2])
                        +dy3*ty1*(u[i][j-1][k][2]
                                 -2.0*u[i][j][k][2]
                                  +u[i][j+1][k][2]);
        rsd[i][j][k][3] += ty3*C3*C4*(flux[i][j+1][k][3]-flux[i][j][k][3])
                        +dy4*ty1*(u[i][j-1][k][3]
                                 -2.0*u[i][j][k][3]
                                  +u[i][j+1][k][3]);
        rsd[i][j][k][4] += ty3*C3*C4*(flux[i][j+1][k][4]-flux[i][j][k][4])
                        +dy5*ty1*(u[i][j-1][k][4]
                                 -2.0*u[i][j][k][4]
                                  +u[i][j+1][k][4]);
      }
      /* Y‐boundary DS */
      for (m = 0; m < 5; m++) {
        rsd[i][1][k][m]   -= dssp*(+5.0*u[i][1][k][m]
                                 -4.0*u[i][2][k][m]
                                  +   u[i][3][k][m]);
        rsd[i][2][k][m]   -= dssp*(-4.0*u[i][1][k][m]
                                 +6.0*u[i][2][k][m]
                                 -4.0*u[i][3][k][m]
                                  +   u[i][4][k][m]);
      }
      /* internal DS in Y */
      jst1 = 3; jend1 = ny-4;
      for (j = jst1; j <= jend1; j++)
        for (m = 0; m < 5; m++)
          rsd[i][j][k][m] -= dssp*(u[i][j-2][k][m]
                                -4.0*u[i][j-1][k][m]
                                +6.0*u[i][j][k][m]
                                -4.0*u[i][j+1][k][m]
                                 +u[i][j+2][k][m]);
      /* Y‐edge DS */
      for (m = 0; m < 5; m++) {
        rsd[i][ny-3][k][m] -= dssp*( u[i][ny-5][k][m]
                                  -4.0*u[i][ny-4][k][m]
                                  +6.0*u[i][ny-3][k][m]
                                  -4.0*u[i][ny-2][k][m]);
        rsd[i][ny-2][k][m] -= dssp*( u[i][ny-4][k][m]
                                  -4.0*u[i][ny-3][k][m]
                                  +5.0*u[i][ny-2][k][m]);
      }
    }
  }

  /*-------------------------------------------------------------------------
    compute fluxes in Z‐direction and add to rsd
  -------------------------------------------------------------------------*/
#pragma omp parallel for collapse(3) private(i,j,k,u41,q) schedule(static)
  for (i = ist; i <= iend; i++) {
    for (j = jst; j <= jend; j++) {
      for (k = 0; k <= nz-1; k++) {
        flux[i][j][k][0] = u[i][j][k][3];
        u41 = u[i][j][k][3] / u[i][j][k][0];
        q = 0.5 * (u[i][j][k][1]*u[i][j][k][1]
                + u[i][j][k][2]*u[i][j][k][2]
                + u[i][j][k][3]*u[i][j][k][3])
            / u[i][j][k][0];
        flux[i][j][k][1] = u[i][j][k][1]*u41;
        flux[i][j][k][2] = u[i][j][k][2]*u41;
        flux[i][j][k][3] = u[i][j][k][3]*u41 + C2*(u[i][j][k][4]-q);
        flux[i][j][k][4] = (C1*u[i][j][k][4] - C2*q)*u41;
      }
    }
  }
#pragma omp parallel for collapse(3) private(i,j,k,m, \
      tmp,u21k,u31k,u41k,u51k, u21km1,u31km1,u41km1,u51km1) schedule(static)
  for (i = ist; i <= iend; i++) {
    for (j = jst; j <= jend; j++) {
      /* flux diff */
      for (k = 1; k <= nz-2; k++) {
        for (m = 0; m < 5; m++) {
          rsd[i][j][k][m] -= tz2*(flux[i][j][k+1][m] - flux[i][j][k-1][m]);
        }
      }
      /* viscous flux and add viscous+DS */
      for (k = 1; k <= nz-1; k++) {
        tmp = 1.0 / u[i][j][k][0];
        u21k = tmp * u[i][j][k][1];
        u31k = tmp * u[i][j][k][2];
        u41k = tmp * u[i][j][k][3];
        u51k = tmp * u[i][j][k][4];
        tmp = 1.0 / u[i][j][k-1][0];
        u21km1 = tmp * u[i][j][k-1][1];
        u31km1 = tmp * u[i][j][k-1][2];
        u41km1 = tmp * u[i][j][k-1][3];
        u51km1 = tmp * u[i][j][k-1][4];
        flux[i][j][k][1] = tz3*(u21k - u21km1);
        flux[i][j][k][2] = tz3*(u31k - u31km1);
        flux[i][j][k][3] = (4.0/3.0)*tz3*(u41k-u41km1);
        flux[i][j][k][4] = 0.5*(1.0-C1*C5)*tz3*
                           ((u21k*u21k+u31k*u31k+u41k*u41k)
                          -(u21km1*u21km1+u31km1*u31km1+u41km1*u41km1))
                          + (1.0/6.0)*tz3*(u41k*u41k-u41km1*u41km1)
                          + C1*C5*tz3*(u51k-u51km1);
      }
      for (k = 1; k <= nz-2; k++) {
        rsd[i][j][k][0] += dz1*tz1*(u[i][j][k-1][0]
                                 -2.0*u[i][j][k][0]
                                  +u[i][j][k+1][0]);
        rsd[i][j][k][1] += tz3*C3*C4*(flux[i][j][k+1][1]-flux[i][j][k][1])
                        +dz2*tz1*(u[i][j][k-1][1]
                                 -2.0*u[i][j][k][1]
                                  +u[i][j][k+1][1]);
        rsd[i][j][k][2] += tz3*C3*C4*(flux[i][j][k+1][2]-flux[i][j][k][2])
                        +dz3*tz1*(u[i][j][k-1][2]
                                 -2.0*u[i][j][k][2]
                                  +u[i][j][k+1][2]);
        rsd[i][j][k][3] += tz3*C3*C4*(flux[i][j][k+1][3]-flux[i][j][k][3])
                        +dz4*tz1*(u[i][j][k-1][3]
                                 -2.0*u[i][j][k][3]
                                  +u[i][j][k+1][3]);
        rsd[i][j][k][4] += tz3*C3*C4*(flux[i][j][k+1][4]-flux[i][j][k][4])
                        +dz5*tz1*(u[i][j][k-1][4]
                                 -2.0*u[i][j][k][4]
                                  +u[i][j][k+1][4]);
      }
      /* Z‐boundary DS (k=1,2) */
      for (m = 0; m < 5; m++) {
        rsd[i][j][1][m]   -= dssp*(+5.0*u[i][j][1][m]
                                 -4.0*u[i][j][2][m]
                                  +   u[i][j][3][m]);
        rsd[i][j][2][m]   -= dssp*(-4.0*u[i][j][1][m]
                                 +6.0*u[i][j][2][m]
                                 -4.0*u[i][j][3][m]
                                  +   u[i][j][4][m]);
      }
      /* Z‐internal DS */
      for (k = 3; k <= nz-4; k++)
        for (m = 0; m < 5; m++)
          rsd[i][j][k][m] -= dssp*(u[i][j][k-2][m]
                                -4.0*u[i][j][k-1][m]
                                +6.0*u[i][j][k][m]
                                -4.0*u[i][j][k+1][m]
                                 +u[i][j][k+2][m]);
      /* Z‐edge DS */
      for (m = 0; m < 5; m++) {
        rsd[i][j][nz-3][m] -= dssp*( u[i][j][nz-5][m]
                                  -4.0*u[i][j][nz-4][m]
                                  +6.0*u[i][j][nz-3][m]
                                  -4.0*u[i][j][nz-2][m]);
        rsd[i][j][nz-2][m] -= dssp*( u[i][j][nz-4][m]
                                  -4.0*u[i][j][nz-3][m]
                                  +5.0*u[i][j][nz-2][m]);
      }
    }
  }
}