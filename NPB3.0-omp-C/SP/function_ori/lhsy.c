static void lhsy(void) {

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c This function computes the left hand side for the three y-factors   
c-------------------------------------------------------------------*/

  double ru1;
  int i, j, k;

/*--------------------------------------------------------------------
c      first fill the lhs for the u-eigenvalue         
c-------------------------------------------------------------------*/
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (k = 1; k <= grid_points[2]-2; k++) {
#pragma omp for  
      for (j = 0; j <= grid_points[1]-1; j++) {
	ru1 = c3c4*rho_i[i][j][k];
	cv[j] = vs[i][j][k];
	rhoq[j] = max(dy3 + con43 * ru1,
		      max(dy5 + c1c5*ru1,
			  max(dymax + ru1,
			      dy1)));
      }
            
#pragma omp for  
      for (j = 1; j <= grid_points[1]-2; j++) {
	lhs[0][i][j][k] =  0.0;
	lhs[1][i][j][k] = -dtty2 * cv[j-1] - dtty1 * rhoq[j-1];
	lhs[2][i][j][k] =  1.0 + c2dtty1 * rhoq[j];
	lhs[3][i][j][k] =  dtty2 * cv[j+1] - dtty1 * rhoq[j+1];
	lhs[4][i][j][k] =  0.0;
      }
    }
  }

/*--------------------------------------------------------------------
c      add fourth order dissipation                             
c-------------------------------------------------------------------*/

  j = 1;
#pragma omp for nowait
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (k = 1; k <= grid_points[2]-2; k++) {

      lhs[2][i][j][k] = lhs[2][i][j][k] + comz5;
      lhs[3][i][j][k] = lhs[3][i][j][k] - comz4;
      lhs[4][i][j][k] = lhs[4][i][j][k] + comz1;
       
      lhs[1][i][j+1][k] = lhs[1][i][j+1][k] - comz4;
      lhs[2][i][j+1][k] = lhs[2][i][j+1][k] + comz6;
      lhs[3][i][j+1][k] = lhs[3][i][j+1][k] - comz4;
      lhs[4][i][j+1][k] = lhs[4][i][j+1][k] + comz1;
    }
  }

#pragma omp for nowait
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 3; j <= grid_points[1]-4; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	lhs[0][i][j][k] = lhs[0][i][j][k] + comz1;
	lhs[1][i][j][k] = lhs[1][i][j][k] - comz4;
	lhs[2][i][j][k] = lhs[2][i][j][k] + comz6;
	lhs[3][i][j][k] = lhs[3][i][j][k] - comz4;
	lhs[4][i][j][k] = lhs[4][i][j][k] + comz1;
      }
    }
  }

  j = grid_points[1]-3;
#pragma omp for  
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (k = 1; k <= grid_points[2]-2; k++) {
      lhs[0][i][j][k] = lhs[0][i][j][k] + comz1;
      lhs[1][i][j][k] = lhs[1][i][j][k] - comz4;
      lhs[2][i][j][k] = lhs[2][i][j][k] + comz6;
      lhs[3][i][j][k] = lhs[3][i][j][k] - comz4;

      lhs[0][i][j+1][k] = lhs[0][i][j+1][k] + comz1;
      lhs[1][i][j+1][k] = lhs[1][i][j+1][k] - comz4;
      lhs[2][i][j+1][k] = lhs[2][i][j+1][k] + comz5;
    }
  }

/*--------------------------------------------------------------------
c      subsequently, do the other two factors                    
c-------------------------------------------------------------------*/
#pragma omp for  
  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
	lhs[0+5][i][j][k]  = lhs[0][i][j][k];
	lhs[1+5][i][j][k]  = lhs[1][i][j][k] - 
	  dtty2 * speed[i][j-1][k];
	lhs[2+5][i][j][k]  = lhs[2][i][j][k];
	lhs[3+5][i][j][k]  = lhs[3][i][j][k] + 
	  dtty2 * speed[i][j+1][k];
	lhs[4+5][i][j][k] = lhs[4][i][j][k];
	lhs[0+10][i][j][k] = lhs[0][i][j][k];
	lhs[1+10][i][j][k] = lhs[1][i][j][k] + 
	  dtty2 * speed[i][j-1][k];
	lhs[2+10][i][j][k] = lhs[2][i][j][k];
	lhs[3+10][i][j][k] = lhs[3][i][j][k] - 
	  dtty2 * speed[i][j+1][k];
	lhs[4+10][i][j][k] = lhs[4][i][j][k];
      }
    }
  }
}