static void pintgr(void) {
  int i, j, k;
  int ibeg, ifin, ifin1;
  int jbeg, jfin, jfin1;
  int iglob, iglob1, iglob2;
  int jglob, jglob1, jglob2;
  double phi1[ISIZ2+2][ISIZ3+2];	
  double phi2[ISIZ2+2][ISIZ3+2];	
  double  frc1, frc2, frc3;
  ibeg = nx;
  ifin = 0;
  iglob1 = -1;
  iglob2 = nx-1;
  if (iglob1 >= ii1 && iglob2 < ii2+nx) ibeg = 0;
  if (iglob1 >= ii1-nx && iglob2 <= ii2) ifin = nx;
  if (ii1 >= iglob1 && ii1 <= iglob2) ibeg = ii1;
  if (ii2 >= iglob1 && ii2 <= iglob2) ifin = ii2;
  jbeg = ny;
  jfin = -1;
  jglob1 = 0;
  jglob2 = ny-1;
  if (jglob1 >= ji1 && jglob2 < ji2+ny) jbeg = 0;
  if (jglob1 > ji1-ny && jglob2 <= ji2) jfin = ny;
  if (ji1 >= jglob1 && ji1 <= jglob2) jbeg = ji1;
  if (ji2 >= jglob1 && ji2 <= jglob2) jfin = ji2;
  ifin1 = ifin;
  jfin1 = jfin;
  if (ifin1 == ii2) ifin1 = ifin -1;
  if (jfin1 == ji2) jfin1 = jfin -1;
  for (i = 0; i <= ISIZ2+1; i++) {
#pragma omp parallel for collapse(2) private(i, k)
    for (k = 0; k <= ISIZ3+1; k++) {
      phi1[i][k] = 0.0;
      phi2[i][k] = 0.0;
    }
  }
  for (i = ibeg; i <= ifin; i++) {
    iglob = i;
#pragma omp parallel for collapse(2) private(i, j, iglob, jglob)
    for (j = jbeg; j <= jfin; j++) {
      jglob = j;
      k = ki1;
      phi1[i][j] = C2*(  u[i][j][k][4]
			- 0.50 * (  pow2(u[i][j][k][1])
				    + pow2(u[i][j][k][2])
				    + pow2(u[i][j][k][3]) )
			/ u[i][j][k][0] );
      k = ki2;
      phi2[i][j] = C2*(  u[i][j][k][4]
			- 0.50 * (  pow2(u[i][j][k][1])
				    + pow2(u[i][j][k][2])
				    + pow2(u[i][j][k][3]) )
			/ u[i][j][k][0] );
    }
  }
  frc1 = 0.0;
  for (i = ibeg; i <= ifin1; i++) {
#pragma omp parallel for collapse(2) private(i, j) reduction(+:frc1)
    for (j = jbeg; j <= jfin1; j++) {
      frc1 = frc1 + (  phi1[i][j]
		       + phi1[i+1][j]
		       + phi1[i][j+1]
		       + phi1[i+1][j+1]
		       + phi2[i][j]
		       + phi2[i+1][j]
		       + phi2[i][j+1]
		       + phi2[i+1][j+1] );
    }
  }
  frc1 = dxi * deta * frc1;
  for (i = 0; i <= ISIZ2+1; i++) {
    for (k = 0; k <= ISIZ3+1; k++) {
#pragma omp parallel for collapse(2) private(i, k)
      phi1[i][k] = 0.0;
      phi2[i][k] = 0.0;
    }
  }
  jglob = jbeg;
  if (jglob == ji1) {
    for (i = ibeg; i <= ifin; i++) {
#pragma omp parallel for collapse(2) private(i, k, iglob)
      iglob = i;
      for (k = ki1; k <= ki2; k++) {
	phi1[i][k] = C2*(  u[i][jbeg][k][4]
			  - 0.50 * (  pow2(u[i][jbeg][k][1])
				      + pow2(u[i][jbeg][k][2])
				      + pow2(u[i][jbeg][k][3]) )
			  / u[i][jbeg][k][0] );
      }
    }
  }
  jglob = jfin;
  if (jglob == ji2) {
    for (i = ibeg; i <= ifin; i++) {
#pragma omp parallel for collapse(2) private(i, k, iglob)
      iglob = i;
      for (k = ki1; k <= ki2; k++) {
	phi2[i][k] = C2*(  u[i][jfin][k][4]
			  - 0.50 * (  pow2(u[i][jfin][k][1])
				      + pow2(u[i][jfin][k][2])
				      + pow2(u[i][jfin][k][3]) )
			  / u[i][jfin][k][0] );
      }
    }
  }
  frc2 = 0.0;
  for (i = ibeg; i <= ifin1; i++) {
#pragma omp parallel for collapse(2) private(i, k) reduction(+:frc2)
    for (k = ki1; k <= ki2-1; k++) {
      frc2 = frc2 + (  phi1[i][k]
		       + phi1[i+1][k]
		       + phi1[i][k+1]
		       + phi1[i+1][k+1]
		       + phi2[i][k]
		       + phi2[i+1][k]
		       + phi2[i][k+1]
		       + phi2[i+1][k+1] );
    }
  }
  frc2 = dxi * dzeta * frc2;
  for (i = 0; i <= ISIZ2+1; i++) {
    for (k = 0; k <= ISIZ3+1; k++) {
#pragma omp parallel for collapse(2) private(i, k)
      phi1[i][k] = 0.0;
      phi2[i][k] = 0.0;
    }
  }
  iglob = ibeg;
  if (iglob == ii1) {
    for (j = jbeg; j <= jfin; j++) {
#pragma omp parallel for collapse(2) private(j, k, jglob)
      jglob = j;
      for (k = ki1; k <= ki2; k++) {
	phi1[j][k] = C2*(  u[ibeg][j][k][4]
			  - 0.50 * (  pow2(u[ibeg][j][k][1])
				      + pow2(u[ibeg][j][k][2])
				      + pow2(u[ibeg][j][k][3]) )
			  / u[ibeg][j][k][0] );
      }
    }
  }
  iglob = ifin;
  if (iglob == ii2) {
    for (j = jbeg; j <= jfin; j++) {
#pragma omp parallel for collapse(2) private(j, k, jglob)
      jglob = j;
      for (k = ki1; k <= ki2; k++) {
	phi2[j][k] = C2*(  u[ifin][j][k][4]
			  - 0.50 * (  pow2(u[ifin][j][k][1])
				      + pow2(u[ifin][j][k][2])
				      + pow2(u[ifin][j][k][3]) )
			  / u[ifin][j][k][0] );
      }
    }
  }
  frc3 = 0.0;
  for (j = jbeg; j <= jfin1; j++) {
#pragma omp parallel for collapse(2) private(j, k) reduction(+:frc3)
    for (k = ki1; k <= ki2-1; k++) {
      frc3 = frc3 + (  phi1[j][k]
		       + phi1[j+1][k]
		       + phi1[j][k+1]
		       + phi1[j+1][k+1]
		       + phi2[j][k]
		       + phi2[j+1][k]
		       + phi2[j][k+1]
		       + phi2[j+1][k+1] );
    }
  }
  frc3 = deta * dzeta * frc3;
  frc = 0.25 * ( frc1 + frc2 + frc3 );
}