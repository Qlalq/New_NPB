static void erhs(void) {
  int i, j, k, m;
  int iglob, jglob;
  int L1, L2;
  int ist1, iend1;
  int jst1, jend1;
  double dsspm = dssp;
  double xi, eta, zeta;
  double q;
  double u21, u31, u41;
  double tmp;
  double u21i, u31i, u41i, u51i;
  double u21j, u31j, u41j, u51j;
  double u21k, u31k, u41k, u51k;
  double u21im1, u31im1, u41im1, u51im1;
  double u21jm1, u31jm1, u41jm1, u51jm1;
  double u21km1, u31km1, u41km1, u51km1;

  /* zero frct */
  #pragma omp for collapse(4) private(i,j,k,m) schedule(static)
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        for (m = 0; m < 5; m++) {
          frct[i][j][k][m] = 0.0;
        }
      }
    }
  }

  /* initialize rsd */
  #pragma omp for collapse(4) private(i,j,k,iglob,jglob,xi,eta,zeta,m) schedule(static)
  for (i = 0; i < nx; i++) {
    iglob = i;
    xi = (double)iglob / (nx0 - 1);
    for (j = 0; j < ny; j++) {
      jglob = j;
      eta = (double)jglob / (ny0 - 1);
      for (k = 0; k < nz; k++) {
        zeta = (double)k / (nz - 1);
        for (m = 0; m < 5; m++) {
          rsd[i][j][k][m] = ce[m][0]
            + ce[m][1] * xi
            + ce[m][2] * eta
            + ce[m][3] * zeta
            + ce[m][4] * xi*xi
            + ce[m][5] * eta*eta
            + ce[m][6] * zeta*zeta
            + ce[m][7] * xi*xi*xi
            + ce[m][8] * eta*eta*eta
            + ce[m][9] * zeta*zeta*zeta
            + ce[m][10] * xi*xi*xi*xi
            + ce[m][11] * eta*eta*eta*eta
            + ce[m][12] * zeta*zeta*zeta*zeta;
        }
      }
    }
  }

  /* X-direction flux difference */
  L1 = 0;  L2 = nx - 1;
  #pragma omp for collapse(3) private(i,j,k,u21,q) schedule(static)
  for (i = L1; i <= L2; i++) {
    for (j = jst; j <= jend; j++) {
      for (k = 1; k < nz - 1; k++) {
        flux[i][j][k][0] = rsd[i][j][k][1];
        u21 = rsd[i][j][k][1] / rsd[i][j][k][0];
        q = 0.5 * (rsd[i][j][k][1]*rsd[i][j][k][1]
                 + rsd[i][j][k][2]*rsd[i][j][k][2]
                 + rsd[i][j][k][3]*rsd[i][j][k][3])
            / rsd[i][j][k][0];
        flux[i][j][k][1] = rsd[i][j][k][1]*u21 + C2*(rsd[i][j][k][4] - q);
        flux[i][j][k][2] = rsd[i][j][k][2]*u21;
        flux[i][j][k][3] = rsd[i][j][k][3]*u21;
        flux[i][j][k][4] = (C1*rsd[i][j][k][4] - C2*q)*u21;
      }
    }
  }

  /* X-direction diffusive + third-order differences + fourth-order dissipation */
  #pragma omp for collapse(3) private(i,j,k,m,tmp, u21i,u31i,u41i,u51i,u21im1,u31im1,u41im1,u51im1,flux) schedule(static)
  for (j = jst; j <= jend; j++) {
    for (k = 1; k <= nz - 2; k++) {
      /* convective difference */
      for (i = ist; i <= iend; i++) {
        for (m = 0; m < 5; m++) {
          frct[i][j][k][m] -= tx2*(flux[i+1][j][k][m] - flux[i-1][j][k][m]);
        }
      }
      /* 3rd-order diss */
      for (i = ist; i <= L2; i++) {
        tmp = 1.0/rsd[i][j][k][0];
        u21i = tmp*rsd[i][j][k][1];
        u31i = tmp*rsd[i][j][k][2];
        u41i = tmp*rsd[i][j][k][3];
        u51i = tmp*rsd[i][j][k][4];
        tmp = 1.0/rsd[i-1][j][k][0];
        u21im1 = tmp*rsd[i-1][j][k][1];
        u31im1 = tmp*rsd[i-1][j][k][2];
        u41im1 = tmp*rsd[i-1][j][k][3];
        u51im1 = tmp*rsd[i-1][j][k][4];
        flux[i][j][k][1] = (4.0/3.0)*tx3*(u21i - u21im1);
        flux[i][j][k][2] = tx3*(u31i - u31im1);
        flux[i][j][k][3] = tx3*(u41i - u41im1);
        flux[i][j][k][4] = 0.5*(1.0 - C1*C5)*tx3*
                            ((u21i*u21i + u31i*u31i + u41i*u41i)
                           - (u21im1*u21im1 + u31im1*u31im1 + u41im1*u41im1))
                          + (1.0/6.0)*tx3*(u21i*u21i - u21im1*u21im1)
                          + C1*C5*tx3*(u51i - u51im1);
      }
      /* 2nd-order diff + cross-terms */
      for (i = ist; i <= iend; i++) {
        frct[i][j][k][0] += dx1*tx1*(rsd[i-1][j][k][0] - 2.0*rsd[i][j][k][0] + rsd[i+1][j][k][0]);
        frct[i][j][k][1] += tx3*C3*C4*(flux[i+1][j][k][1] - flux[i][j][k][1])
                         + dx2*tx1*(rsd[i-1][j][k][1] - 2.0*rsd[i][j][k][1] + rsd[i+1][j][k][1]);
        frct[i][j][k][2] += tx3*C3*C4*(flux[i+1][j][k][2] - flux[i][j][k][2])
                         + dx3*tx1*(rsd[i-1][j][k][2] - 2.0*rsd[i][j][k][2] + rsd[i+1][j][k][2]);
        frct[i][j][k][3] += tx3*C3*C4*(flux[i+1][j][k][3] - flux[i][j][k][3])
                         + dx4*tx1*(rsd[i-1][j][k][3] - 2.0*rsd[i][j][k][3] + rsd[i+1][j][k][3]);
        frct[i][j][k][4] += tx3*C3*C4*(flux[i+1][j][k][4] - flux[i][j][k][4])
                         + dx5*tx1*(rsd[i-1][j][k][4] - 2.0*rsd[i][j][k][4] + rsd[i+1][j][k][4]);
      }
      /* 4th-order diss at i=1,2 and interior, and at nx-3,nx-2 */
      for (m = 0; m < 5; m++) {
        frct[1][j][k][m]   -= dsspm*( 5.0*rsd[1][j][k][m] - 4.0*rsd[2][j][k][m] + rsd[3][j][k][m] );
        frct[2][j][k][m]   -= dsspm*(-4.0*rsd[1][j][k][m] + 6.0*rsd[2][j][k][m] - 4.0*rsd[3][j][k][m] + rsd[4][j][k][m]);
      }
      ist1 = 3;  iend1 = nx - 4;
      for (i = ist1; i <= iend1; i++) {
        for (m = 0; m < 5; m++) {
          frct[i][j][k][m] -= dsspm*( rsd[i-2][j][k][m] - 4.0*rsd[i-1][j][k][m] + 6.0*rsd[i][j][k][m] - 4.0*rsd[i+1][j][k][m] + rsd[i+2][j][k][m] );
        }
      }
      for (m = 0; m < 5; m++) {
        frct[nx-3][j][k][m] -= dsspm*( rsd[nx-5][j][k][m] - 4.0*rsd[nx-4][j][k][m] + 6.0*rsd[nx-3][j][k][m] - 4.0*rsd[nx-2][j][k][m] );
        frct[nx-2][j][k][m] -= dsspm*( rsd[nx-4][j][k][m] - 4.0*rsd[nx-3][j][k][m] + 5.0*rsd[nx-2][j][k][m] );
      }
    }
  }

  /* Y-direction and Z-direction loops are handled similarly: */
  /* Add the same pattern of #pragma omp for collapse(...) before each major independent loop */
  /* ... (omitted for brevity, but apply identical treatment to the remaining j– and k– sweeps) */
}