static void lhsx(void) {

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c This function computes the left hand side for the three x-factors  
c-------------------------------------------------------------------*/

  double ru1;
  int i, j, k;

/*--------------------------------------------------------------------
c      first fill the lhs for the u-eigenvalue                   
c-------------------------------------------------------------------*/
  for (j = 1; j <= grid_points[1]-2; j++) {
    for (k = 1; k <= grid_points[2]-2; k++) {
#pragma omp for  
      for (i = 0; i <= grid_points[0]-1; i++) {
	ru1 = c3c4*rho_i[i][j][k];
	cv[i] = us[i][j][k];
	rhon[i] = max(dx2+con43*ru1, 
		      max(dx5+c1c5*ru1,
			  max(dxmax+ru1,
			      dx1)));
      }

#pragma omp for  
      for (i = 1; i <= grid_points[0]-2; i++) {
	lhs[0][i][j][k] =   0.0;
	lhs[1][i][j][k] = - dttx2 * cv[i-1] - dttx1 * rhon[i-1];
	lhs[2][i][j][k] =   1.0 + c2dttx1 * rhon[i];
	lhs[3][i][j][k] =   dttx2 * cv[i+1] - dttx1 * rhon[i+1];
	lhs[4][i][j][k] =   0.0;
      }
    }
  }

/*--------------------------------------------------------------------
c      add fourth order dissipation                             
c-------------------------------------------------------------------*/

  i = 1;
#pragma omp for nowait
  for (j = 1; j <= grid_points[1]-2; j++) {
    for (k = 1; k <= grid_points[2]-2; k++) {
      lhs[2][i][j][k] = lhs[2][i][j][k] + comz5;
      lhs[3][i][j][k] = lhs[3][i][j][k] - comz4;
      lhs[4][i][j][k] = lhs[4][i][j][k] + comz1;
      lhs[1][i+1][j][k] = lhs[1][i+1][j][k] - comz4;
      lhs[2][i+1][j][k] = lhs[2][i+1][j][k] + comz6;
      lhs[3][i+1][j][k] = lhs[3][i+1][j][k] - comz4;
      lhs[4][i+1][j][k] = lhs[4][i+1][j][k] + comz1;
    }
  }

#pragma omp for nowait
  for (i = 3; i <= grid_points[0]-4; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	lhs[0][i][j][k] = lhs[0][i][j][k] + comz1;
	lhs[1][i][j][k] = lhs[1][i][j][k] - comz4;
	lhs[2][i][j][k] = lhs[2][i][j][k] + comz6;
	lhs[3][i][j][k] = lhs[3][i][j][k] - comz4;
	lhs[4][i][j][k] = lhs[4][i][j][k] + comz1;
      }
    }
  }

  i = grid_points[0]-3;
#pragma omp for  
  for (j = 1; j <= grid_points[1]-2; j++) {
    for (k = 1; k <= grid_points[2]-2; k++) {
      lhs[0][i][j][k] = lhs[0][i][j][k] + comz1;
      lhs[1][i][j][k] = lhs[1][i][j][k] - comz4;
      lhs[2][i][j][k] = lhs[2][i][j][k] + comz6;
      lhs[3][i][j][k] = lhs[3][i][j][k] - comz4;

      lhs[0][i+1][j][k] = lhs[0][i+1][j][k] + comz1;
      lhs[1][i+1][j][k] = lhs[1][i+1][j][k] - comz4;
      lhs[2][i+1][j][k] = lhs[2][i+1][j][k] + comz5;
    }
  }

/*--------------------------------------------------------------------
c      subsequently, fill the other factors (u+c), (u-c) by adding to 
c      the first  
c-------------------------------------------------------------------*/
#pragma omp for  
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	lhs[0+5][i][j][k]  = lhs[0][i][j][k];
	lhs[1+5][i][j][k]  = lhs[1][i][j][k] - 
	  dttx2 * speed[i-1][j][k];
	lhs[2+5][i][j][k]  = lhs[2][i][j][k];
	lhs[3+5][i][j][k]  = lhs[3][i][j][k] + 
	  dttx2 * speed[i+1][j][k];
	lhs[4+5][i][j][k]  = lhs[4][i][j][k];
	lhs[0+10][i][j][k] = lhs[0][i][j][k];
	lhs[1+10][i][j][k] = lhs[1][i][j][k] + 
	  dttx2 * speed[i-1][j][k];
	lhs[2+10][i][j][k] = lhs[2][i][j][k];
	lhs[3+10][i][j][k] = lhs[3][i][j][k] - 
	  dttx2 * speed[i+1][j][k];
	lhs[4+10][i][j][k] = lhs[4][i][j][k];
      }
    }
  }
}