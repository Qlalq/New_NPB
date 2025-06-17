static void exact_rhs(void) {

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c compute the right hand side based on exact solution
c-------------------------------------------------------------------*/

  double dtemp[5], xi, eta, zeta, dtpp;
  int m, i, j, k, ip1, im1, jp1, jm1, km1, kp1;

/*--------------------------------------------------------------------
c      initialize                                  
c-------------------------------------------------------------------*/
  for (m = 0; m < 5; m++) {
    for (i = 0; i <= grid_points[0]-1; i++) {
      for (j = 0; j <= grid_points[1]-1; j++) {
	for (k= 0; k <= grid_points[2]-1; k++) {
	  forcing[m][i][j][k] = 0.0;
	}
      }
    }
  }

/*--------------------------------------------------------------------
c      xi-direction flux differences                      
c-------------------------------------------------------------------*/
  for (k = 1; k <= grid_points[2]-2; k++) {
    zeta = (double)k * dnzm1;
    for (j = 1; j <= grid_points[1]-2; j++) {
      eta = (double)j * dnym1;

      for (i = 0; i <= grid_points[0]-1; i++) {
	xi = (double)i * dnxm1;

	exact_solution(xi, eta, zeta, dtemp);
	for (m = 0; m < 5; m++) {
	  ue[m][i] = dtemp[m];
	}

	dtpp = 1.0 / dtemp[0];

	for (m = 1; m < 5; m++) {
	  buf[m][i] = dtpp * dtemp[m];
	}

	cuf[i] = buf[1][i] * buf[1][i];
	buf[0][i] = cuf[i] + buf[2][i] * buf[2][i] + buf[3][i] * buf[3][i];
	q[i] = 0.5 * (buf[1][i]*ue[1][i] + buf[2][i]*ue[2][i]
		      + buf[3][i]*ue[3][i]);
      }
 
      for (i = 1; i <= grid_points[0]-2; i++) {
	im1 = i-1;
	ip1 = i+1;

	forcing[0][i][j][k] = forcing[0][i][j][k] -
	  tx2*( ue[1][ip1]-ue[1][im1] )+
	  dx1tx1*(ue[0][ip1]-2.0*ue[0][i]+ue[0][im1]);

	forcing[1][i][j][k] = forcing[1][i][j][k]
	  - tx2 * ((ue[1][ip1]*buf[1][ip1]+c2*(ue[4][ip1]-q[ip1]))-
                   (ue[1][im1]*buf[1][im1]+c2*(ue[4][im1]-q[im1])))+
	  xxcon1*(buf[1][ip1]-2.0*buf[1][i]+buf[1][im1])+
	  dx2tx1*( ue[1][ip1]-2.0* ue[1][i]+ue[1][im1]);

	forcing[2][i][j][k] = forcing[2][i][j][k]
	  - tx2 * (ue[2][ip1]*buf[1][ip1]-ue[2][im1]*buf[1][im1])+
	  xxcon2*(buf[2][ip1]-2.0*buf[2][i]+buf[2][im1])+
	  dx3tx1*( ue[2][ip1]-2.0*ue[2][i] +ue[2][im1]);
                  
	forcing[3][i][j][k] = forcing[3][i][j][k]
	  - tx2*(ue[3][ip1]*buf[1][ip1]-ue[3][im1]*buf[1][im1])+
	  xxcon2*(buf[3][ip1]-2.0*buf[3][i]+buf[3][im1])+
	  dx4tx1*( ue[3][ip1]-2.0* ue[3][i]+ ue[3][im1]);

	forcing[4][i][j][k] = forcing[4][i][j][k]
	  - tx2*(buf[1][ip1]*(c1*ue[4][ip1]-c2*q[ip1])-
		 buf[1][im1]*(c1*ue[4][im1]-c2*q[im1]))+
	  0.5*xxcon3*(buf[0][ip1]-2.0*buf[0][i]+
		      buf[0][im1])+
	  xxcon4*(cuf[ip1]-2.0*cuf[i]+cuf[im1])+
	  xxcon5*(buf[4][ip1]-2.0*buf[4][i]+buf[4][im1])+
	  dx5tx1*( ue[4][ip1]-2.0* ue[4][i]+ ue[4][im1]);
      }

/*--------------------------------------------------------------------
c            Fourth-order dissipation                         
c-------------------------------------------------------------------*/
      for (m = 0; m < 5; m++) {
	i = 1;
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (5.0*ue[m][i] - 4.0*ue[m][i+1] +ue[m][i+2]);
	i = 2;
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (-4.0*ue[m][i-1] + 6.0*ue[m][i] -
 	    4.0*ue[m][i+1] +     ue[m][i+2]);
      }

      for (m = 0; m < 5; m++) {
	for (i = 3; i <= grid_points[0]-4; i++) {
	  forcing[m][i][j][k] = forcing[m][i][j][k] - dssp*
	    (ue[m][i-2] - 4.0*ue[m][i-1] +
	     6.0*ue[m][i] - 4.0*ue[m][i+1] + ue[m][i+2]);
	}
      }

      for (m = 0; m < 5; m++) {
	i = grid_points[0]-3;
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (ue[m][i-2] - 4.0*ue[m][i-1] +
	   6.0*ue[m][i] - 4.0*ue[m][i+1]);
	i = grid_points[0]-2;
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (ue[m][i-2] - 4.0*ue[m][i-1] + 5.0*ue[m][i]);
      }
    }
  }

/*--------------------------------------------------------------------
c  eta-direction flux differences             
c-------------------------------------------------------------------*/
  for (k = 1; k <= grid_points[2]-2; k++) {
    zeta = (double)k * dnzm1;
    for (i = 1; i <= grid_points[0]-2; i++) {
      xi = (double)i * dnxm1;

      for (j = 0; j <= grid_points[1]-1; j++) {
	eta = (double)j * dnym1;

	exact_solution(xi, eta, zeta, dtemp);
	for (m = 0; m < 5; m++) {
	  ue[m][j] = dtemp[m];
	}
	dtpp = 1.0/dtemp[0];

	for (m = 1; m < 5; m++) {
	  buf[m][j] = dtpp * dtemp[m];
	}

	cuf[j]   = buf[2][j] * buf[2][j];
	buf[0][j] = cuf[j] + buf[1][j] * buf[1][j] + 
	  buf[3][j] * buf[3][j];
	q[j] = 0.5*(buf[1][j]*ue[1][j] + buf[2][j]*ue[2][j] +
		    buf[3][j]*ue[3][j]);
      }

      for (j = 1; j <= grid_points[1]-2; j++) {
	jm1 = j-1;
	jp1 = j+1;
                  
	forcing[0][i][j][k] = forcing[0][i][j][k] -
	  ty2*( ue[2][jp1]-ue[2][jm1] )+
	  dy1ty1*(ue[0][jp1]-2.0*ue[0][j]+ue[0][jm1]);

	forcing[1][i][j][k] = forcing[1][i][j][k]
	  - ty2*(ue[1][jp1]*buf[2][jp1]-ue[1][jm1]*buf[2][jm1])+
	  yycon2*(buf[1][jp1]-2.0*buf[1][j]+buf[1][jm1])+
	  dy2ty1*( ue[1][jp1]-2.0* ue[1][j]+ ue[1][jm1]);

	forcing[2][i][j][k] = forcing[2][i][j][k]
	  - ty2*((ue[2][jp1]*buf[2][jp1]+c2*(ue[4][jp1]-q[jp1]))-
		 (ue[2][jm1]*buf[2][jm1]+c2*(ue[4][jm1]-q[jm1])))+
	  yycon1*(buf[2][jp1]-2.0*buf[2][j]+buf[2][jm1])+
	  dy3ty1*( ue[2][jp1]-2.0*ue[2][j] +ue[2][jm1]);

	forcing[3][i][j][k] = forcing[3][i][j][k]
	  - ty2*(ue[3][jp1]*buf[2][jp1]-ue[3][jm1]*buf[2][jm1])+
	  yycon2*(buf[3][jp1]-2.0*buf[3][j]+buf[3][jm1])+
	  dy4ty1*( ue[3][jp1]-2.0*ue[3][j]+ ue[3][jm1]);

	forcing[4][i][j][k] = forcing[4][i][j][k]
	  - ty2*(buf[2][jp1]*(c1*ue[4][jp1]-c2*q[jp1])-
		 buf[2][jm1]*(c1*ue[4][jm1]-c2*q[jm1]))+
	  0.5*yycon3*(buf[0][jp1]-2.0*buf[0][j]+
		      buf[0][jm1])+
	  yycon4*(cuf[jp1]-2.0*cuf[j]+cuf[jm1])+
	  yycon5*(buf[4][jp1]-2.0*buf[4][j]+buf[4][jm1])+
	  dy5ty1*(ue[4][jp1]-2.0*ue[4][j]+ue[4][jm1]);
      }

/*--------------------------------------------------------------------
c            Fourth-order dissipation                      
c-------------------------------------------------------------------*/
      for (m = 0; m < 5; m++) {
	j = 1;
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (5.0*ue[m][j] - 4.0*ue[m][j+1] +ue[m][j+2]);
	j = 2;
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (-4.0*ue[m][j-1] + 6.0*ue[m][j] -
	   4.0*ue[m][j+1] +       ue[m][j+2]);
      }

      for (m = 0; m < 5; m++) {
	for (j = 3; j <= grid_points[1]-4; j++) {
	  forcing[m][i][j][k] = forcing[m][i][j][k] - dssp*
	    (ue[m][j-2] - 4.0*ue[m][j-1] +
	     6.0*ue[m][j] - 4.0*ue[m][j+1] + ue[m][j+2]);
	}
      }

      for (m = 0; m < 5; m++) {
	j = grid_points[1]-3;
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (ue[m][j-2] - 4.0*ue[m][j-1] +
	   6.0*ue[m][j] - 4.0*ue[m][j+1]);
	j = grid_points[1]-2;
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (ue[m][j-2] - 4.0*ue[m][j-1] + 5.0*ue[m][j]);

      }
    }
  }

/*--------------------------------------------------------------------
c      zeta-direction flux differences                      
c-------------------------------------------------------------------*/
  for (j = 1; j <= grid_points[1]-2; j++) {
    eta = (double)j * dnym1;
    for (i = 1; i <= grid_points[0]-2; i++) {
      xi = (double)i * dnxm1;

      for (k = 0; k <= grid_points[2]-1; k++) {
	zeta = (double)k * dnzm1;

	exact_solution(xi, eta, zeta, dtemp);
	for (m = 0; m < 5; m++) {
	  ue[m][k] = dtemp[m];
	}

	dtpp = 1.0/dtemp[0];

	for (m = 1; m < 5; m++) {
	  buf[m][k] = dtpp * dtemp[m];
	}

	cuf[k] = buf[3][k] * buf[3][k];
	buf[0][k] = cuf[k] + buf[1][k] * buf[1][k] + 
	  buf[2][k] * buf[2][k];
	q[k] = 0.5*(buf[1][k]*ue[1][k] + buf[2][k]*ue[2][k] +
		    buf[3][k]*ue[3][k]);
      }

      for (k = 1; k <= grid_points[2]-2; k++) {
	km1 = k-1;
	kp1 = k+1;
                  
	forcing[0][i][j][k] = forcing[0][i][j][k] -
	  tz2*( ue[3][kp1]-ue[3][km1] )+
	  dz1tz1*(ue[0][kp1]-2.0*ue[0][k]+ue[0][km1]);

	forcing[1][i][j][k] = forcing[1][i][j][k]
	  - tz2 * (ue[1][kp1]*buf[3][kp1]-ue[1][km1]*buf[3][km1])+
	  zzcon2*(buf[1][kp1]-2.0*buf[1][k]+buf[1][km1])+
	  dz2tz1*( ue[1][kp1]-2.0* ue[1][k]+ ue[1][km1]);

	forcing[2][i][j][k] = forcing[2][i][j][k]
	  - tz2 * (ue[2][kp1]*buf[3][kp1]-ue[2][km1]*buf[3][km1])+
	  zzcon2*(buf[2][kp1]-2.0*buf[2][k]+buf[2][km1])+
	  dz3tz1*(ue[2][kp1]-2.0*ue[2][k]+ue[2][km1]);

	forcing[3][i][j][k] = forcing[3][i][j][k]
	  - tz2 * ((ue[3][kp1]*buf[3][kp1]+c2*(ue[4][kp1]-q[kp1]))-
		   (ue[3][km1]*buf[3][km1]+c2*(ue[4][km1]-q[km1])))+
	  zzcon1*(buf[3][kp1]-2.0*buf[3][k]+buf[3][km1])+
	  dz4tz1*( ue[3][kp1]-2.0*ue[3][k] +ue[3][km1]);

	forcing[4][i][j][k] = forcing[4][i][j][k]
	  - tz2 * (buf[3][kp1]*(c1*ue[4][kp1]-c2*q[kp1])-
		   buf[3][km1]*(c1*ue[4][km1]-c2*q[km1]))+
	  0.5*zzcon3*(buf[0][kp1]-2.0*buf[0][k]
		      +buf[0][km1])+
	  zzcon4*(cuf[kp1]-2.0*cuf[k]+cuf[km1])+
	  zzcon5*(buf[4][kp1]-2.0*buf[4][k]+buf[4][km1])+
	  dz5tz1*( ue[4][kp1]-2.0*ue[4][k]+ ue[4][km1]);
      }

/*--------------------------------------------------------------------
c            Fourth-order dissipation                        
c-------------------------------------------------------------------*/
      for (m = 0; m < 5; m++) {
	k = 1;
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (5.0*ue[m][k] - 4.0*ue[m][k+1] +ue[m][k+2]);
	k = 2;
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (-4.0*ue[m][k-1] + 6.0*ue[m][k] -
	   4.0*ue[m][k+1] +       ue[m][k+2]);
      }

      for (m = 0; m < 5; m++) {
	for (k = 3; k <= grid_points[2]-4; k++) {
	  forcing[m][i][j][k] = forcing[m][i][j][k] - dssp*
	    (ue[m][k-2] - 4.0*ue[m][k-1] +
	     6.0*ue[m][k] - 4.0*ue[m][k+1] + ue[m][k+2]);
	}
      }

      for (m = 0; m < 5; m++) {
	k = grid_points[2]-3;
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (ue[m][k-2] - 4.0*ue[m][k-1] +
	   6.0*ue[m][k] - 4.0*ue[m][k+1]);
	k = grid_points[2]-2;
	forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
	  (ue[m][k-2] - 4.0*ue[m][k-1] + 5.0*ue[m][k]);
      }
    }
  }

/*--------------------------------------------------------------------
c now change the sign of the forcing function, 
c-------------------------------------------------------------------*/
  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = 1; k <= grid_points[2]-2; k++) {
	  forcing[m][i][j][k] = -1.0 * forcing[m][i][j][k];
	}
      }
    }
  }
}