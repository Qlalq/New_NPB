static void compute_rhs(void) {
  int i, j, k, m;
  double aux, rho_inv, uijk, up1, um1, vijk, vp1, vm1,
    wijk, wp1, wm1;

  //---------------------------------------------------------------------
  // Step 1: Calculate intermediate values based on u
  // These calculations are independent for each (i, j, k) point.
  // This loop nest is fully parallelizable over i, j, and k.
  // No data dependencies between iterations.
  //---------------------------------------------------------------------
  for (i = 0; i <= grid_points[0]-1; i++) {
    for (j = 0; j <= grid_points[1]-1; j++) {
      for (k = 0; k <= grid_points[2]-1; k++) {
	rho_inv = 1.0/u[0][i][j][k];
	rho_i[i][j][k] = rho_inv;
	us[i][j][k] = u[1][i][j][k] * rho_inv;
	vs[i][j][k] = u[2][i][j][k] * rho_inv;
	ws[i][j][k] = u[3][i][j][k] * rho_inv;
	square[i][j][k] = 0.5* (u[1][i][j][k]*u[1][i][j][k] +
				u[2][i][j][k]*u[2][i][j][k] +
				u[3][i][j][k]*u[3][i][j][k] ) * rho_inv;
	qs[i][j][k] = square[i][j][k] * rho_inv;
	aux = c1c2*rho_inv* (u[4][i][j][k] - square[i][j][k]);
	aux = sqrt(aux); // Assuming aux >= 0
	speed[i][j][k] = aux;
	ainv[i][j][k]  = 1.0/aux; // Assuming aux != 0
      }
    }
  }

  //---------------------------------------------------------------------
  // Step 2: Initialize rhs with forcing term
  // These assignments are independent for each (m, i, j, k) point.
  // This loop nest is fully parallelizable over m, i, j, and k.
  // No data dependencies.
  //---------------------------------------------------------------------
  for (m = 0; m < 5; m++) {
    for (i = 0; i <= grid_points[0]-1; i++) {
      for (j = 0; j <= grid_points[1]-1; j++) {
	for (k = 0; k <= grid_points[2]-1; k++) {
	  rhs[m][i][j][k] = forcing[m][i][j][k];
	}
      }
    }
  }

  //---------------------------------------------------------------------
  // Step 3: Add x-direction contributions to rhs
  // This block accumulates contributions from x-directional central differences
  // and artificial viscosity stencils. Each loop nest within this block
  // updates a distinct set of rhs[m][i][j][k] points based on spatial neighbors.
  // The loops within each nest are parallelizable over spatial indices
  // (i, j, k) and 'm'. This block must be completed before the y-direction
  // contributions are added to avoid race conditions on rhs.
  //---------------------------------------------------------------------

  // Add 3-point stencil central difference contributions in x-direction (interior points)
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	uijk = us[i][j][k];
	up1  = us[i+1][j][k];
	um1  = us[i-1][j][k];
	rhs[0][i][j][k] = rhs[0][i][j][k] + dx1tx1 *
	  (u[0][i+1][j][k] - 2.0*u[0][i][j][k] +
	   u[0][i-1][j][k]) -
	  tx2 * (u[1][i+1][j][k] - u[1][i-1][j][k]);
	rhs[1][i][j][k] = rhs[1][i][j][k] + dx2tx1 *
	  (u[1][i+1][j][k] - 2.0*u[1][i][j][k] +
	   u[1][i-1][j][k]) +
	  xxcon2*con43 * (up1 - 2.0*uijk + um1) -
	  tx2 * (u[1][i+1][j][k]*up1 -
		 u[1][i-1][j][k]*um1 +
		 (u[4][i+1][j][k]- square[i+1][j][k]-
		  u[4][i-1][j][k]+ square[i-1][j][k])*
		 c2);
	rhs[2][i][j][k] = rhs[2][i][j][k] + dx3tx1 *
	  (u[2][i+1][j][k] - 2.0*u[2][i][j][k] +
	   u[2][i-1][j][k]) +
	  xxcon2 * (vs[i+1][j][k] - 2.0*vs[i][j][k] +
		    vs[i-1][j][k]) -
	  tx2 * (u[2][i+1][j][k]*up1 -
		 u[2][i-1][j][k]*um1);
	rhs[3][i][j][k] = rhs[3][i][j][k] + dx4tx1 *
	  (u[3][i+1][j][k] - 2.0*u[3][i][j][k] +
	   u[3][i-1][j][k]) +
	  xxcon2 * (ws[i+1][j][k] - 2.0*ws[i][j][k] +
		    ws[i-1][j][k]) -
	  tx2 * (u[3][i+1][j][k]*up1 -
		 u[3][i-1][j][k]*um1);
	rhs[4][i][j][k] = rhs[4][i][j][k] + dx5tx1 *
	  (u[4][i+1][j][k] - 2.0*u[4][i][j][k] +
	   u[4][i-1][j][k]) +
	  xxcon3 * (qs[i+1][j][k] - 2.0*qs[i][j][k] +
		    qs[i-1][j][k]) +
	  xxcon4 * (up1*up1 -       2.0*uijk*uijk +
		    um1*um1) +
	  xxcon5 * (u[4][i+1][j][k]*rho_i[i+1][j][k] -
		    2.0*u[4][i][j][k]*rho_i[i][j][k] +
		    u[4][i-1][j][k]*rho_i[i-1][j][k]) -
	  tx2 * ( (c1*u[4][i+1][j][k] -
		   c2*square[i+1][j][k])*up1 -
		  (c1*u[4][i-1][j][k] -
		   c2*square[i-1][j][k])*um1 );
      }
    }
  }

  // Add 5-point stencil dssp contributions in x-direction (boundary and interior points)
  // i = 1 boundary
  i = 1;
  for (m = 0; m < 5; m++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	rhs[m][i][j][k] = rhs[m][i][j][k]- dssp *
	  ( 5.0*u[m][i][j][k] - 4.0*u[m][i+1][j][k] +
	    u[m][i+2][j][k]);
      }
    }
  }
  // i = 2 boundary
  i = 2;
  for (m = 0; m < 5; m++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
	  (-4.0*u[m][i-1][j][k] + 6.0*u[m][i][j][k] -
	   4.0*u[m][i+1][j][k] + u[m][i+2][j][k]);
      }
    }
  }
  // i = 3 to grid_points[0]-4 interior
  for (m = 0; m < 5; m++) {
    for (i = 3*1; i <= grid_points[0]-3*1-1; i++) { // Assuming 3*1 is just 3 and similar for upper bound
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = 1; k <= grid_points[2]-2; k++) {
	  rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
	    (  u[m][i-2][j][k] - 4.0*u[m][i-1][j][k] +
	       6.0*u[m][i][j][k] - 4.0*u[m][i+1][j][k] +
	       u[m][i+2][j][k] );
	}
      }
    }
  }
  // i = grid_points[0]-3 boundary
  i = grid_points[0]-3;
  for (m = 0; m < 5; m++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
	  ( u[m][i-2][j][k] - 4.0*u[m][i-1][j][k] +
	    6.0*u[m][i][j][k] - 4.0*u[m][i+1][j][k] );
      }
    }
  }
  // i = grid_points[0]-2 boundary
  i = grid_points[0]-2;
  for (m = 0; m < 5; m++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
	  ( u[m][i-2][j][k] - 4.0*u[m][i-1][j][k] +
	    5.0*u[m][i][j][k] );
      }
    }
  }

  //---------------------------------------------------------------------
  // Step 4: Add y-direction contributions to rhs
  // This block accumulates contributions from y-directional central differences
  // and artificial viscosity stencils. Similar parallelization strategy as
  // the x-direction block. This block must be completed before the z-direction
  // contributions are added.
  //---------------------------------------------------------------------

  // Add 3-point stencil central difference contributions in y-direction (interior points)
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	vijk = vs[i][j][k];
	vp1  = vs[i][j+1][k];
	vm1  = vs[i][j-1][k];
	rhs[0][i][j][k] = rhs[0][i][j][k] + dy1ty1 *
	  (u[0][i][j+1][k] - 2.0*u[0][i][j][k] +
	   u[0][i][j-1][k]) -
	  ty2 * (u[2][i][j+1][k] - u[2][i][j-1][k]);
	rhs[1][i][j][k] = rhs[1][i][j][k] + dy2ty1 *
	  (u[1][i][j+1][k] - 2.0*u[1][i][j][k] +
	   u[1][i][j-1][k]) +
	  yycon2 * (us[i][j+1][k] - 2.0*us[i][j][k] +
		    us[i][j-1][k]) -
	  ty2 * (u[1][i][j+1][k]*vp1 -
		 u[1][i][j-1][k]*vm1);
	rhs[2][i][j][k] = rhs[2][i][j][k] + dy3ty1 *
	  (u[2][i][j+1][k] - 2.0*u[2][i][j][k] +
	   u[2][i][j-1][k]) +
	  yycon2*con43 * (vp1 - 2.0*vijk + vm1) -
	  ty2 * (u[2][i][j+1][k]*vp1 -
		 u[2][i][j-1][k]*vm1 +
		 (u[4][i][j+1][k] - square[i][j+1][k] -
		  u[4][i][j-1][k] + square[i][j-1][k])
		 *c2);
	rhs[3][i][j][k] = rhs[3][i][j][k] + dy4ty1 *
	  (u[3][i][j+1][k] - 2.0*u[3][i][j][k] +
	   u[3][i][j-1][k]) +
	  yycon2 * (ws[i][j+1][k] - 2.0*ws[i][j][k] +
		    ws[i][j-1][k]) -
	  ty2 * (u[3][i][j+1][k]*vp1 -
		 u[3][i][j-1][k]*vm1);
	rhs[4][i][j][k] = rhs[4][i][j][k] + dy5ty1 *
	  (u[4][i][j+1][k] - 2.0*u[4][i][j][k] +
	   u[4][i][j-1][k]) +
	  yycon3 * (qs[i][j+1][k] - 2.0*qs[i][j][k] +
		    qs[i][j-1][k]) +
	  yycon4 * (vp1*vp1       - 2.0*vijk*vijk +
		    vm1*vm1) +
	  yycon5 * (u[4][i][j+1][k]*rho_i[i][j+1][k] -
		    2.0*u[4][i][j][k]*rho_i[i][j][k] +
		    u[4][i][j-1][k]*rho_i[i][j-1][k]) -
	  ty2 * ((c1*u[4][i][j+1][k] -
		  c2*square[i][j+1][k]) * vp1 -
		 (c1*u[4][i][j-1][k] -
		  c2*square[i][j-1][k]) * vm1);
      }
    }
  }

  // Add 5-point stencil dssp contributions in y-direction (boundary and interior points)
  // j = 1 boundary
  j = 1;
  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	rhs[m][i][j][k] = rhs[m][i][j][k]- dssp *
	  ( 5.0*u[m][i][j][k] - 4.0*u[m][i][j+1][k] +
	    u[m][i][j+2][k]);
      }
    }
  }
  // j = 2 boundary
  j = 2;
  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
	  (-4.0*u[m][i][j-1][k] + 6.0*u[m][i][j][k] -
	   4.0*u[m][i][j+1][k] + u[m][i][j+2][k]);
      }
    }
  }
  // j = 3 to grid_points[1]-4 interior
  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 3*1; j <= grid_points[1]-3*1-1; j++) { // Assuming 3*1 is just 3 and similar for upper bound
	for (k = 1; k <= grid_points[2]-2; k++) {
	  rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
	    (  u[m][i][j-2][k] - 4.0*u[m][i][j-1][k] +
	       6.0*u[m][i][j][k] - 4.0*u[m][i][j+1][k] +
	       u[m][i][j+2][k] );
	}
      }
    }
  }
  // j = grid_points[1]-3 boundary
  j = grid_points[1]-3;
  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
	  ( u[m][i][j-2][k] - 4.0*u[m][i][j-1][k] +
	    6.0*u[m][i][j][k] - 4.0*u[m][i][j+1][k] );
      }
    }
  }
  // j = grid_points[1]-2 boundary
  j = grid_points[1]-2;
  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
	  ( u[m][i][j-2][k] - 4.0*u[m][i][j-1][k] +
	    5.0*u[m][i][j][k] );
      }
    }
  }

  //---------------------------------------------------------------------
  // Step 5: Add z-direction contributions to rhs
  // This block accumulates contributions from z-directional central differences
  // and artificial viscosity stencils. Similar parallelization strategy as
  // the x and y-direction blocks.
  //---------------------------------------------------------------------

  // Add 3-point stencil central difference contributions in z-direction (interior points)
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	wijk = ws[i][j][k];
	wp1  = ws[i][j][k+1];
	wm1  = ws[i][j][k-1];
	rhs[0][i][j][k] = rhs[0][i][j][k] + dz1tz1 *
	  (u[0][i][j][k+1] - 2.0*u[0][i][j][k] +
	   u[0][i][j][k-1]) -
	  tz2 * (u[3][i][j][k+1] - u[3][i][j][k-1]);
	rhs[1][i][j][k] = rhs[1][i][j][k] + dz2tz1 *
	  (u[1][i][j][k+1] - 2.0*u[1][i][j][k] +
	   u[1][i][j][k-1]) +
	  zzcon2 * (us[i][j][k+1] - 2.0*us[i][j][k] +
		    us[i][j][k-1]) -
	  tz2 * (u[1][i][j][k+1]*wp1 -
		 u[1][i][j][k-1]*wm1);
	rhs[2][i][j][k] = rhs[2][i][j][k] + dz3tz1 *
	  (u[2][i][j][k+1] - 2.0*u[2][i][j][k] +
	   u[2][i][j][k-1]) +
	  zzcon2 * (vs[i][j][k+1] - 2.0*vs[i][j][k] +
		    vs[i][j][k-1]) -
	  tz2 * (u[2][i][j][k+1]*wp1 -
		 u[2][i][j][k-1]*wm1);
	rhs[3][i][j][k] = rhs[3][i][j][k] + dz4tz1 *
	  (u[3][i][j][k+1] - 2.0*u[3][i][j][k] +
	   u[3][i][j][k-1]) +
	  zzcon2*con43 * (wp1 - 2.0*wijk + wm1) -
	  tz2 * (u[3][i][j][k+1]*wp1 -
		 u[3][i][j][k-1]*wm1 +
		 (u[4][i][j][k+1] - square[i][j][k+1] -
		  u[4][i][j][k-1] + square[i][j][k-1])
		 *c2);
	rhs[4][i][j][k] = rhs[4][i][j][k] + dz5tz1 *
	  (u[4][i][j][k+1] - 2.0*u[4][i][j][k] +
	   u[4][i][j][k-1]) +
	  zzcon3 * (qs[i][j][k+1] - 2.0*qs[i][j][k] +
		    qs[i][j][k-1]) +
	  zzcon4 * (wp1*wp1 - 2.0*wijk*wijk +
		    wm1*wm1) +
	  zzcon5 * (u[4][i][j][k+1]*rho_i[i][j][k+1] -
		    2.0*u[4][i][j][k]*rho_i[i][j][k] +
		    u[4][i][j][k-1]*rho_i[i][j][k-1]) -
	  tz2 * ( (c1*u[4][i][j][k+1] -
		   c2*square[i][j][k+1])*wp1 -
		  (c1*u[4][i][j][k-1] -
		   c2*square[i][j][k-1])*wm1);
      }
    }
  }

  // Add 5-point stencil dssp contributions in z-direction (boundary and interior points)
  // k = 1 boundary
  k = 1;
  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	rhs[m][i][j][k] = rhs[m][i][j][k]- dssp *
	  ( 5.0*u[m][i][j][k] - 4.0*u[m][i][j][k+1] +
	    u[m][i][j][k+2]);
      }
    }
  }
  // k = 2 boundary
  k = 2;
  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = 1; k <= grid_points[2]-2; k++) { // Note: original loop included k from 1..size-2, which seems off for a k=2 boundary, but preserving original code.
	rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
	  (-4.0*u[m][i][j][k-1] + 6.0*u[m][i][j][k] -
	   4.0*u[m][i][j][k+1] + u[m][i][j][k+2]);
      }
    }
  }
  // k = 3 to grid_points[2]-4 interior
  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = 3*1; k <= grid_points[2]-3*1-1; k++) { // Assuming 3*1 is just 3 and similar for upper bound
	  rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
	    (  u[m][i][j][k-2] - 4.0*u[m][i][j][k-1] +
	       6.0*u[m][i][j][k] - 4.0*u[m][i][j][k+1] +
	       u[m][i][j][k+2] );
	}
      }
    }
  }
  // k = grid_points[2]-3 boundary
  k = grid_points[2]-3;
  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = 1; k <= grid_points[2]-2; k++) { // Note: original loop included k from 1..size-2, which seems off for a k=grid_points[2]-3 boundary, but preserving original code.
	  rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
	    ( u[m][i][j][k-2] - 4.0*u[m][i][j][k-1] +
	      6.0*u[m][i][j][k] - 4.0*u[m][i][j][k+1] );
	}
      }
    }
  }
  // k = grid_points[2]-2 boundary
  k = grid_points[2]-2;
  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = 1; k <= grid_points[2]-2; k++) { // Note: original loop included k from 1..size-2, which seems off for a k=grid_points[2]-2 boundary, but preserving original code.
	  rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
	    ( u[m][i][j][k-2] - 4.0*u[m][i][j][k-1] +
	      5.0*u[m][i][j][k] );
	}
      }
    }
  }

  //---------------------------------------------------------------------
  // Step 6: Apply time step factor (dt)
  // This scaling operation is independent for each (m, i, j, k) point
  // where rhs was modified in the previous steps.
  // This loop is fully parallelizable over m, i, j, and k.
  // No data dependencies.
  // Note: Original range is 1..size-2 for i,j,k, consistent with stencil regions.
  //---------------------------------------------------------------------
  for (m = 0; m < 5; m++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
	for (k = 1; k <= grid_points[2]-2; k++) {
	  rhs[m][i][j][k] = rhs[m][i][j][k] * dt;
	}
      }
    }
  }
}