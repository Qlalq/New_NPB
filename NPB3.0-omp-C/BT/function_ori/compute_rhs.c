static void compute_rhs(void) {

  int i, j, k, m;
  double rho_inv, uijk, up1, um1, vijk, vp1, vm1, wijk, wp1, wm1;

/*--------------------------------------------------------------------
c     compute the reciprocal of density, and the kinetic energy, 
c     and the speed of sound.
c-------------------------------------------------------------------*/
#pragma omp for nowait
  for (i = 0; i < grid_points[0]; i++) {
    for (j = 0; j < grid_points[1]; j++) {
      for (k = 0; k < grid_points[2]; k++) {
	rho_inv = 1.0/u[i][j][k][0];
	rho_i[i][j][k] = rho_inv;
	us[i][j][k] = u[i][j][k][1] * rho_inv;
	vs[i][j][k] = u[i][j][k][2] * rho_inv;
	ws[i][j][k] = u[i][j][k][3] * rho_inv;
	square[i][j][k] = 0.5 * (u[i][j][k][1]*u[i][j][k][1] + 
				 u[i][j][k][2]*u[i][j][k][2] +
				 u[i][j][k][3]*u[i][j][k][3] ) * rho_inv;
	qs[i][j][k] = square[i][j][k] * rho_inv;
      }
    }
  }

/*--------------------------------------------------------------------
c copy the exact forcing term to the right hand side;  because 
c this forcing term is known, we can store it on the whole grid
c including the boundary                   
c-------------------------------------------------------------------*/

#pragma omp for
  for (i = 0; i < grid_points[0]; i++) {
    for (j = 0; j < grid_points[1]; j++) {
      for (k = 0; k < grid_points[2]; k++) {
	for (m = 0; m < 5; m++) {
	  rhs[i][j][k][m] = forcing[i][j][k][m];
	}
      }
    }
  }

/*--------------------------------------------------------------------
c     compute xi-direction fluxes 
c-------------------------------------------------------------------*/
#pragma omp for
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	uijk = us[i][j][k];
	up1  = us[i+1][j][k];
	um1  = us[i-1][j][k];

	rhs[i][j][k][0] = rhs[i][j][k][0] + dx1tx1 * 
	  (u[i+1][j][k][0] - 2.0*u[i][j][k][0] + 
	   u[i-1][j][k][0]) -
	  tx2 * (u[i+1][j][k][1] - u[i-1][j][k][1]);

	rhs[i][j][k][1] = rhs[i][j][k][1] + dx2tx1 * 
	  (u[i+1][j][k][1] - 2.0*u[i][j][k][1] + 
	   u[i-1][j][k][1]) +
	  xxcon2*con43 * (up1 - 2.0*uijk + um1) -
	  tx2 * (u[i+1][j][k][1]*up1 - 
		 u[i-1][j][k][1]*um1 +
		 (u[i+1][j][k][4]- square[i+1][j][k]-
		  u[i-1][j][k][4]+ square[i-1][j][k])*
		 c2);

	rhs[i][j][k][2] = rhs[i][j][k][2] + dx3tx1 * 
	  (u[i+1][j][k][2] - 2.0*u[i][j][k][2] +
	   u[i-1][j][k][2]) +
	  xxcon2 * (vs[i+1][j][k] - 2.0*vs[i][j][k] +
		    vs[i-1][j][k]) -
	  tx2 * (u[i+1][j][k][2]*up1 - 
		 u[i-1][j][k][2]*um1);

	rhs[i][j][k][3] = rhs[i][j][k][3] + dx4tx1 * 
	  (u[i+1][j][k][3] - 2.0*u[i][j][k][3] +
	   u[i-1][j][k][3]) +
	  xxcon2 * (ws[i+1][j][k] - 2.0*ws[i][j][k] +
		    ws[i-1][j][k]) -
	  tx2 * (u[i+1][j][k][3]*up1 - 
		 u[i-1][j][k][3]*um1);

	rhs[i][j][k][4] = rhs[i][j][k][4] + dx5tx1 * 
	  (u[i+1][j][k][4] - 2.0*u[i][j][k][4] +
	   u[i-1][j][k][4]) +
	  xxcon3 * (qs[i+1][j][k] - 2.0*qs[i][j][k] +
		    qs[i-1][j][k]) +
	  xxcon4 * (up1*up1 -       2.0*uijk*uijk + 
		    um1*um1) +
	  xxcon5 * (u[i+1][j][k][4]*rho_i[i+1][j][k] - 
		    2.0*u[i][j][k][4]*rho_i[i][j][k] +
		    u[i-1][j][k][4]*rho_i[i-1][j][k]) -
	  tx2 * ( (c1*u[i+1][j][k][4] - 
		   c2*square[i+1][j][k])*up1 -
		  (c1*u[i-1][j][k][4] - 
		   c2*square[i-1][j][k])*um1 );
      }
    }
  }

/*--------------------------------------------------------------------
c     add fourth order xi-direction dissipation               
c-------------------------------------------------------------------*/
  i = 1;
#pragma omp for nowait
  for (j = 1; j < grid_points[1]-1; j++) {
    for (k = 1; k < grid_points[2]-1; k++) {
      for (m = 0; m < 5; m++) {
	rhs[i][j][k][m] = rhs[i][j][k][m]- dssp * 
	  ( 5.0*u[i][j][k][m] - 4.0*u[i+1][j][k][m] +
	    u[i+2][j][k][m]);
      }
    }
  }

  i = 2;
#pragma omp for nowait
  for (j = 1; j < grid_points[1]-1; j++) {
    for (k = 1; k < grid_points[2]-1; k++) {
      for (m = 0; m < 5; m++) {
	rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * 
	  (-4.0*u[i-1][j][k][m] + 6.0*u[i][j][k][m] -
	   4.0*u[i+1][j][k][m] + u[i+2][j][k][m]);
      }
    }
  }

#pragma omp for nowait
  for (i = 3; i < grid_points[0]-3; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	for (m = 0; m < 5; m++) {
	  rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * 
	    (  u[i-2][j][k][m] - 4.0*u[i-1][j][k][m] + 
	       6.0*u[i][j][k][m] - 4.0*u[i+1][j][k][m] + 
	       u[i+2][j][k][m] );
	}
      }
    }
  }
         
  i = grid_points[0]-3;
#pragma omp for nowait
  for (j = 1; j < grid_points[1]-1; j++) {
    for (k = 1; k < grid_points[2]-1; k++) {
      for (m = 0; m < 5; m++) {
	rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
	  ( u[i-2][j][k][m] - 4.0*u[i-1][j][k][m] + 
	    6.0*u[i][j][k][m] - 4.0*u[i+1][j][k][m] );
      }
    }
  }

  i = grid_points[0]-2;
#pragma omp for
  for (j = 1; j < grid_points[1]-1; j++) {
    for (k = 1; k < grid_points[2]-1; k++) {
      for (m = 0; m < 5; m++) {
	rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
	  ( u[i-2][j][k][m] - 4.*u[i-1][j][k][m] +
	    5.0*u[i][j][k][m] );
      }
    }
  }

/*--------------------------------------------------------------------
c     compute eta-direction fluxes 
c-------------------------------------------------------------------*/
#pragma omp for
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	vijk = vs[i][j][k];
	vp1  = vs[i][j+1][k];
	vm1  = vs[i][j-1][k];
	rhs[i][j][k][0] = rhs[i][j][k][0] + dy1ty1 * 
	  (u[i][j+1][k][0] - 2.0*u[i][j][k][0] + 
	   u[i][j-1][k][0]) -
	  ty2 * (u[i][j+1][k][2] - u[i][j-1][k][2]);
	rhs[i][j][k][1] = rhs[i][j][k][1] + dy2ty1 * 
	  (u[i][j+1][k][1] - 2.0*u[i][j][k][1] + 
	   u[i][j-1][k][1]) +
	  yycon2 * (us[i][j+1][k] - 2.0*us[i][j][k] + 
		    us[i][j-1][k]) -
	  ty2 * (u[i][j+1][k][1]*vp1 - 
		 u[i][j-1][k][1]*vm1);
	rhs[i][j][k][2] = rhs[i][j][k][2] + dy3ty1 * 
	  (u[i][j+1][k][2] - 2.0*u[i][j][k][2] + 
	   u[i][j-1][k][2]) +
	  yycon2*con43 * (vp1 - 2.0*vijk + vm1) -
	  ty2 * (u[i][j+1][k][2]*vp1 - 
		 u[i][j-1][k][2]*vm1 +
		 (u[i][j+1][k][4] - square[i][j+1][k] - 
		  u[i][j-1][k][4] + square[i][j-1][k])
		 *c2);
	rhs[i][j][k][3] = rhs[i][j][k][3] + dy4ty1 * 
	  (u[i][j+1][k][3] - 2.0*u[i][j][k][3] + 
	   u[i][j-1][k][3]) +
	  yycon2 * (ws[i][j+1][k] - 2.0*ws[i][j][k] + 
		    ws[i][j-1][k]) -
	  ty2 * (u[i][j+1][k][3]*vp1 - 
		 u[i][j-1][k][3]*vm1);
	rhs[i][j][k][4] = rhs[i][j][k][4] + dy5ty1 * 
	  (u[i][j+1][k][4] - 2.0*u[i][j][k][4] + 
	   u[i][j-1][k][4]) +
	  yycon3 * (qs[i][j+1][k] - 2.0*qs[i][j][k] + 
		    qs[i][j-1][k]) +
	  yycon4 * (vp1*vp1       - 2.0*vijk*vijk + 
		    vm1*vm1) +
	  yycon5 * (u[i][j+1][k][4]*rho_i[i][j+1][k] - 
		    2.0*u[i][j][k][4]*rho_i[i][j][k] +
		    u[i][j-1][k][4]*rho_i[i][j-1][k]) -
	  ty2 * ((c1*u[i][j+1][k][4] - 
		  c2*square[i][j+1][k]) * vp1 -
		 (c1*u[i][j-1][k][4] - 
		  c2*square[i][j-1][k]) * vm1);
      }
    }
  }

/*--------------------------------------------------------------------
c     add fourth order eta-direction dissipation         
c-------------------------------------------------------------------*/
  j = 1;
#pragma omp for nowait
  for (i = 1; i < grid_points[0]-1; i++) {
    for (k = 1; k < grid_points[2]-1; k++) {
      for (m = 0; m < 5; m++) {
	rhs[i][j][k][m] = rhs[i][j][k][m]- dssp * 
	  ( 5.0*u[i][j][k][m] - 4.0*u[i][j+1][k][m] +
	    u[i][j+2][k][m]);
      }
    }
  }

  j = 2;
#pragma omp for nowait
  for (i = 1; i < grid_points[0]-1; i++) {
    for (k = 1; k < grid_points[2]-1; k++) {
      for (m = 0; m < 5; m++) {
	rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * 
	  (-4.0*u[i][j-1][k][m] + 6.0*u[i][j][k][m] -
	   4.0*u[i][j+1][k][m] + u[i][j+2][k][m]);
      }
    }
  }

#pragma omp for nowait
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 3; j < grid_points[1]-3; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	for (m = 0; m < 5; m++) {
	  rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * 
	    (  u[i][j-2][k][m] - 4.0*u[i][j-1][k][m] + 
	       6.0*u[i][j][k][m] - 4.0*u[i][j+1][k][m] + 
	       u[i][j+2][k][m] );
	}
      }
    }
  }
         
  j = grid_points[1]-3;
#pragma omp for nowait
  for (i = 1; i < grid_points[0]-1; i++) {
    for (k = 1; k < grid_points[2]-1; k++) {
      for (m = 0; m < 5; m++) {
	rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
	  ( u[i][j-2][k][m] - 4.0*u[i][j-1][k][m] + 
	    6.0*u[i][j][k][m] - 4.0*u[i][j+1][k][m] );
      }
    }
  }

  j = grid_points[1]-2;
#pragma omp for
  for (i = 1; i < grid_points[0]-1; i++) {
    for (k = 1; k < grid_points[2]-1; k++) {
      for (m = 0; m < 5; m++) {
	rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
	  ( u[i][j-2][k][m] - 4.*u[i][j-1][k][m] +
	    5.*u[i][j][k][m] );
      }
    }
  }

/*--------------------------------------------------------------------
c     compute zeta-direction fluxes 
c-------------------------------------------------------------------*/
#pragma omp for
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	wijk = ws[i][j][k];
	wp1  = ws[i][j][k+1];
	wm1  = ws[i][j][k-1];

	rhs[i][j][k][0] = rhs[i][j][k][0] + dz1tz1 * 
	  (u[i][j][k+1][0] - 2.0*u[i][j][k][0] + 
	   u[i][j][k-1][0]) -
	  tz2 * (u[i][j][k+1][3] - u[i][j][k-1][3]);
	rhs[i][j][k][1] = rhs[i][j][k][1] + dz2tz1 * 
	  (u[i][j][k+1][1] - 2.0*u[i][j][k][1] + 
	   u[i][j][k-1][1]) +
	  zzcon2 * (us[i][j][k+1] - 2.0*us[i][j][k] + 
		    us[i][j][k-1]) -
	  tz2 * (u[i][j][k+1][1]*wp1 - 
		 u[i][j][k-1][1]*wm1);
	rhs[i][j][k][2] = rhs[i][j][k][2] + dz3tz1 * 
	  (u[i][j][k+1][2] - 2.0*u[i][j][k][2] + 
	   u[i][j][k-1][2]) +
	  zzcon2 * (vs[i][j][k+1] - 2.0*vs[i][j][k] + 
		    vs[i][j][k-1]) -
	  tz2 * (u[i][j][k+1][2]*wp1 - 
		 u[i][j][k-1][2]*wm1);
	rhs[i][j][k][3] = rhs[i][j][k][3] + dz4tz1 * 
	  (u[i][j][k+1][3] - 2.0*u[i][j][k][3] + 
	   u[i][j][k-1][3]) +
	  zzcon2*con43 * (wp1 - 2.0*wijk + wm1) -
	  tz2 * (u[i][j][k+1][3]*wp1 - 
		 u[i][j][k-1][3]*wm1 +
		 (u[i][j][k+1][4] - square[i][j][k+1] - 
		  u[i][j][k-1][4] + square[i][j][k-1])
		 *c2);
	rhs[i][j][k][4] = rhs[i][j][k][4] + dz5tz1 * 
	  (u[i][j][k+1][4] - 2.0*u[i][j][k][4] + 
	   u[i][j][k-1][4]) +
	  zzcon3 * (qs[i][j][k+1] - 2.0*qs[i][j][k] + 
		    qs[i][j][k-1]) +
	  zzcon4 * (wp1*wp1 - 2.0*wijk*wijk + 
		    wm1*wm1) +
	  zzcon5 * (u[i][j][k+1][4]*rho_i[i][j][k+1] - 
		    2.0*u[i][j][k][4]*rho_i[i][j][k] +
		    u[i][j][k-1][4]*rho_i[i][j][k-1]) -
	  tz2 * ( (c1*u[i][j][k+1][4] - 
		   c2*square[i][j][k+1])*wp1 -
		  (c1*u[i][j][k-1][4] - 
		   c2*square[i][j][k-1])*wm1);
      }
    }
  }

/*--------------------------------------------------------------------
c     add fourth order zeta-direction dissipation                
c-------------------------------------------------------------------*/
  k = 1;
#pragma omp for nowait
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (m = 0; m < 5; m++) {
	rhs[i][j][k][m] = rhs[i][j][k][m]- dssp * 
	  ( 5.0*u[i][j][k][m] - 4.0*u[i][j][k+1][m] +
	    u[i][j][k+2][m]);
      }
    }
  }

  k = 2;
#pragma omp for nowait
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (m = 0; m < 5; m++) {
	rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * 
	  (-4.0*u[i][j][k-1][m] + 6.0*u[i][j][k][m] -
	   4.0*u[i][j][k+1][m] + u[i][j][k+2][m]);
      }
    }
  }

#pragma omp for nowait
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 3; k < grid_points[2]-3; k++) {
	for (m = 0; m < 5; m++) {
	  rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * 
	    (  u[i][j][k-2][m] - 4.0*u[i][j][k-1][m] + 
	       6.0*u[i][j][k][m] - 4.0*u[i][j][k+1][m] + 
	       u[i][j][k+2][m] );
	}
      }
    }
  }
         
  k = grid_points[2]-3;
#pragma omp for nowait
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (m = 0; m < 5; m++) {
	rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
	  ( u[i][j][k-2][m] - 4.0*u[i][j][k-1][m] + 
	    6.0*u[i][j][k][m] - 4.0*u[i][j][k+1][m] );
      }
    }
  }

  k = grid_points[2]-2;
#pragma omp for
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (m = 0; m < 5; m++) {
	rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
	  ( u[i][j][k-2][m] - 4.0*u[i][j][k-1][m] +
	    5.0*u[i][j][k][m] );
      }
    }
  }

#pragma omp for
  for (j = 1; j < grid_points[1]-1; j++) {
    for (k = 1; k < grid_points[2]-1; k++) {
      for (m = 0; m < 5; m++) {
	for (i = 1; i < grid_points[0]-1; i++) {
	  rhs[i][j][k][m] = rhs[i][j][k][m] * dt;
	}
      }
    }
  }
}