static void initialize(void) {

#pragma omp parallel
{
/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c     This subroutine initializes the field variable u using 
c     tri-linear transfinite interpolation of the boundary values     
c-------------------------------------------------------------------*/

  int i, j, k, m, ix, iy, iz;
  double xi, eta, zeta, Pface[2][3][5], Pxi, Peta, Pzeta, temp[5];

/*--------------------------------------------------------------------
c  Later (in compute_rhs) we compute 1/u for every element. A few of 
c  the corner elements are not used, but it convenient (and faster) 
c  to compute the whole thing with a simple loop. Make sure those 
c  values are nonzero by initializing the whole thing here. 
c-------------------------------------------------------------------*/
#pragma omp for
  for (i = 0; i < IMAX; i++) {
    for (j = 0; j < IMAX; j++) {
      for (k = 0; k < IMAX; k++) {
	for (m = 0; m < 5; m++) {
	  u[i][j][k][m] = 1.0;
	}
      }
    }
  }

/*--------------------------------------------------------------------
c     first store the "interpolated" values everywhere on the grid    
c-------------------------------------------------------------------*/

#pragma omp for
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      
      for (k = 0; k < grid_points[2]; k++) {
	zeta = (double)k * dnzm1;
                  
	for (ix = 0; ix < 2; ix++) {
	  exact_solution((double)ix, eta, zeta, 
                         &(Pface[ix][0][0]));
	}

	for (iy = 0; iy < 2; iy++) {
	  exact_solution(xi, (double)iy , zeta, 
                         &Pface[iy][1][0]);
	}

	for (iz = 0; iz < 2; iz++) {
	  exact_solution(xi, eta, (double)iz,   
                         &Pface[iz][2][0]);
	}

	for (m = 0; m < 5; m++) {
	  Pxi   = xi   * Pface[1][0][m] + 
	    (1.0-xi)   * Pface[0][0][m];
	  Peta  = eta  * Pface[1][1][m] + 
	    (1.0-eta)  * Pface[0][1][m];
	  Pzeta = zeta * Pface[1][2][m] + 
	    (1.0-zeta) * Pface[0][2][m];
                     
	  u[i][j][k][m] = Pxi + Peta + Pzeta - 
	    Pxi*Peta - Pxi*Pzeta - Peta*Pzeta + 
	    Pxi*Peta*Pzeta;
	}
      }
    }
  }

/*--------------------------------------------------------------------
c     now store the exact values on the boundaries        
c-------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c     west face                                                  
c-------------------------------------------------------------------*/
  i = 0;
  xi = 0.0;
#pragma omp for nowait
  for (j = 0; j < grid_points[1]; j++) {
    eta = (double)j * dnym1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m];
      }
    }
  }

/*--------------------------------------------------------------------
c     east face                                                      
c-------------------------------------------------------------------*/

  i = grid_points[0]-1;
  xi = 1.0;
#pragma omp for
  for (j = 0; j < grid_points[1]; j++) {
    eta = (double)j * dnym1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m];
      }
    }
  }

/*--------------------------------------------------------------------
c     south face                                                 
c-------------------------------------------------------------------*/
  j = 0;
  eta = 0.0;
#pragma omp for nowait
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m];
      }
    }
  }

/*--------------------------------------------------------------------
c     north face                                    
c-------------------------------------------------------------------*/
  j = grid_points[1]-1;
  eta = 1.0;
#pragma omp for
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m];
      }
    }
  }

/*--------------------------------------------------------------------
c     bottom face                                       
c-------------------------------------------------------------------*/
  k = 0;
  zeta = 0.0;
#pragma omp for nowait
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i *dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m];
      }
    }
  }

/*--------------------------------------------------------------------
c     top face     
c-------------------------------------------------------------------*/
  k = grid_points[2]-1;
  zeta = 1.0;
#pragma omp for
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m];
      }
    }
  }
}

}