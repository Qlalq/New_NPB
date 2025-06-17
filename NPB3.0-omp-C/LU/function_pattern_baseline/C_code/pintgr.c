static void pintgr(void) {
  int i, j, k;
  int ibeg, ifin, ifin1;
  int jbeg, jfin, jfin1;
  int iglob, iglob1, iglob2;
  int jglob, jglob1, jglob2;
  double phi1[ISIZ2+2][ISIZ3+2];
  double phi2[ISIZ2+2][ISIZ3+2];
  double frc1 = 0.0, frc2 = 0.0, frc3 = 0.0;

  /* compute index bounds */
  ibeg = nx;  ifin = 0;
  iglob1 = -1;  iglob2 = nx-1;
  if (iglob1 >= ii1 && iglob2 < ii2+nx) ibeg = 0;
  if (iglob1 >= ii1-nx && iglob2 <= ii2) ifin = nx;
  if (ii1 >= iglob1 && ii1 <= iglob2) ibeg = ii1;
  if (ii2 >= iglob1 && ii2 <= iglob2) ifin = ii2;

  jbeg = ny;  jfin = -1;
  jglob1 = 0;  jglob2 = ny-1;
  if (jglob1 >= ji1 && jglob2 < ji2+ny) jbeg = 0;
  if (jglob1 > ji1-ny && jglob2 <= ji2) jfin = ny;
  if (ji1 >= jglob1 && ji1 <= jglob2) jbeg = ji1;
  if (ji2 >= jglob1 && ji2 <= jglob2) jfin = ji2;

  ifin1 = ifin;  jfin1 = jfin;
  if (ifin1 == ii2) ifin1 = ifin - 1;
  if (jfin1 == ji2) jfin1 = jfin - 1;

#pragma omp parallel default(shared) private(i,j,k,iglob,jglob) reduction(+:frc1,frc2,frc3)
  {
    /* zero phi arrays */
    #pragma omp for collapse(2) schedule(static)
    for (i = 0; i <= ISIZ2+1; i++) {
      for (k = 0; k <= ISIZ3+1; k++) {
        phi1[i][k] = 0.0;
        phi2[i][k] = 0.0;
      }
    }

    /* compute phi1, phi2 on i–j plane */
    #pragma omp for collapse(2) schedule(static)
    for (i = ibeg; i <= ifin; i++) {
      iglob = i;
      for (j = jbeg; j <= jfin; j++) {
        jglob = j;
        k = ki1;
        phi1[i][j] = C2*(  u[i][j][k][4]
                          - 0.50 * ( pow2(u[i][j][k][1])
                                     + pow2(u[i][j][k][2])
                                     + pow2(u[i][j][k][3]) )
                          / u[i][j][k][0] );
        k = ki2;
        phi2[i][j] = C2*(  u[i][j][k][4]
                          - 0.50 * ( pow2(u[i][j][k][1])
                                     + pow2(u[i][j][k][2])
                                     + pow2(u[i][j][k][3]) )
                          / u[i][j][k][0] );
      }
    }

    /* integrate over i–j to get frc1 */
    #pragma omp for collapse(2) schedule(static)
    for (i = ibeg; i <= ifin1; i++) {
      for (j = jbeg; j <= jfin1; j++) {
        frc1 += (  phi1[i][j]
                  + phi1[i+1][j]
                  + phi1[i][j+1]
                  + phi1[i+1][j+1]
                  + phi2[i][j]
                  + phi2[i+1][j]
                  + phi2[i][j+1]
                  + phi2[i+1][j+1] );
      }
    }

    /* reset phi */
    #pragma omp for collapse(2) schedule(static)
    for (i = 0; i <= ISIZ2+1; i++) {
      for (k = 0; k <= ISIZ3+1; k++) {
        phi1[i][k] = 0.0;
        phi2[i][k] = 0.0;
      }
    }

    /* compute phi1 along j = jbeg boundary */
    #pragma omp for schedule(static)
    for (i = ibeg; i <= ifin; i++) {
      if (jbeg == ji1) {
        iglob = i;
        for (k = ki1; k <= ki2; k++) {
          phi1[i][k] = C2*(  u[i][jbeg][k][4]
                            - 0.50 * ( pow2(u[i][jbeg][k][1])
                                       + pow2(u[i][jbeg][k][2])
                                       + pow2(u[i][jbeg][k][3]) )
                            / u[i][jbeg][k][0] );
        }
      }
    }

    /* compute phi2 along j = jfin boundary */
    #pragma omp for schedule(static)
    for (i = ibeg; i <= ifin; i++) {
      if (jfin == ji2) {
        iglob = i;
        for (k = ki1; k <= ki2; k++) {
          phi2[i][k] = C2*(  u[i][jfin][k][4]
                            - 0.50 * ( pow2(u[i][jfin][k][1])
                                       + pow2(u[i][jfin][k][2])
                                       + pow2(u[i][jfin][k][3]) )
                            / u[i][jfin][k][0] );
        }
      }
    }

    /* integrate over i–k to get frc2 */
    #pragma omp for collapse(2) schedule(static)
    for (i = ibeg; i <= ifin1; i++) {
      for (k = ki1; k <= ki2-1; k++) {
        frc2 += (  phi1[i][k]
                  + phi1[i+1][k]
                  + phi1[i][k+1]
                  + phi1[i+1][k+1]
                  + phi2[i][k]
                  + phi2[i+1][k]
                  + phi2[i][k+1]
                  + phi2[i+1][k+1] );
      }
    }

    /* reset phi again */
    #pragma omp for collapse(2) schedule(static)
    for (i = 0; i <= ISIZ2+1; i++) {
      for (k = 0; k <= ISIZ3+1; k++) {
        phi1[i][k] = 0.0;
        phi2[i][k] = 0.0;
      }
    }

    /* compute phi1 along i = ibeg boundary */
    #pragma omp for schedule(static)
    for (j = jbeg; j <= jfin; j++) {
      if (ibeg == ii1) {
        jglob = j;
        for (k = ki1; k <= ki2; k++) {
          phi1[j][k] = C2*(  u[ibeg][j][k][4]
                            - 0.50 * ( pow2(u[ibeg][j][k][1])
                                       + pow2(u[ibeg][j][k][2])
                                       + pow2(u[ibeg][j][k][3]) )
                            / u[ibeg][j][k][0] );
        }
      }
    }

    /* compute phi2 along i = ifin boundary */
    #pragma omp for schedule(static)
    for (j = jbeg; j <= jfin; j++) {
      if (ifin == ii2) {
        jglob = j;
        for (k = ki1; k <= ki2; k++) {
          phi2[j][k] = C2*(  u[ifin][j][k][4]
                            - 0.50 * ( pow2(u[ifin][j][k][1])
                                       + pow2(u[ifin][j][k][2])
                                       + pow2(u[ifin][j][k][3]) )
                            / u[ifin][j][k][0] );
        }
      }
    }

    /* integrate over j–k to get frc3 */
    #pragma omp for collapse(2) schedule(static)
    for (j = jbeg; j <= jfin1; j++) {
      for (k = ki1; k <= ki2-1; k++) {
        frc3 += (  phi1[j][k]
                  + phi1[j+1][k]
                  + phi1[j][k+1]
                  + phi1[j+1][k+1]
                  + phi2[j][k]
                  + phi2[j+1][k]
                  + phi2[j][k+1]
                  + phi2[j+1][k+1] );
      }
    }
  } /* end omp parallel */

  frc1 *= dxi * deta;
  frc2 *= dxi * dzeta;
  frc3 *= deta * dzeta;
  frc = 0.25 * (frc1 + frc2 + frc3);
}