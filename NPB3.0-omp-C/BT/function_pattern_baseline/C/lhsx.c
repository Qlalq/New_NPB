static void lhsx(void) {
  int i, j, k;
#pragma omp parallel for collapse(2) private(i, tmp1, tmp2, tmp3)
  for (j = 1; j < grid_points[1]-1; j++) {
    for (k = 1; k < grid_points[2]-1; k++) {
      for (i = 0; i < grid_points[0]; i++) {
	tmp1 = 1.0 / u[i][j][k][0];
	tmp2 = tmp1 * tmp1;
	tmp3 = tmp1 * tmp2;
	fjac[ i][ j][ k][0][0] = 0.0;
	fjac[ i][ j][ k][0][1] = 1.0;
	fjac[ i][ j][ k][0][2] = 0.0;
	fjac[ i][ j][ k][0][3] = 0.0;
	fjac[ i][ j][ k][0][4] = 0.0;
	fjac[ i][ j][ k][1][0] = -(u[i][j][k][1] * tmp2 * 
				    u[i][j][k][1])
	  + c2 * 0.50 * (u[i][j][k][1] * u[i][j][k][1]
		       + u[i][j][k][2] * u[i][j][k][2]
		       + u[i][j][k][3] * u[i][j][k][3] ) * tmp2;
	fjac[i][j][k][1][1] = ( 2.0 - c2 )
	  * ( u[i][j][k][1] / u[i][j][k][0] );
	fjac[i][j][k][1][2] = - c2 * ( u[i][j][k][2] * tmp1 );
	fjac[i][j][k][1][3] = - c2 * ( u[i][j][k][3] * tmp1 );
	fjac[i][j][k][1][4] = c2;
	fjac[i][j][k][2][0] = - ( u[i][j][k][1]*u[i][j][k][2] ) * tmp2;
	fjac[i][j][k][2][1] = u[i][j][k][2] * tmp1;
	fjac[i][j][k][2][2] = u[i][j][k][1] * tmp1;
	fjac[i][j][k][2][3] = 0.0;
	fjac[i][j][k][2][4] = 0.0;
	fjac[i][j][k][3][0] = - ( u[i][j][k][1]*u[i][j][k][3] ) * tmp2;
	fjac[i][j][k][3][1] = u[i][j][k][3] * tmp1;
	fjac[i][j][k][3][2] = 0.0;
	fjac[i][j][k][3][3] = u[i][j][k][1] * tmp1;
	fjac[i][j][k][3][4] = 0.0;
	fjac[i][j][k][4][0] = ( c2 * ( u[i][j][k][1] * u[i][j][k][1]
				     + u[i][j][k][2] * u[i][j][k][2]
				     + u[i][j][k][3] * u[i][j][k][3] ) * tmp2
				- c1 * ( u[i][j][k][4] * tmp1 ) )
	  * ( u[i][j][k][1] * tmp1 );
	fjac[i][j][k][4][1] = c1 *  u[i][j][k][4] * tmp1 
	  - 0.50 * c2
	  * (  3.0*u[i][j][k][1]*u[i][j][k][1]
	       + u[i][j][k][2]*u[i][j][k][2]
	       + u[i][j][k][3]*u[i][j][k][3] ) * tmp2;
	fjac[i][j][k][4][2] = - c2 * ( u[i][j][k][2]*u[i][j][k][1] )
	  * tmp2;
	fjac[i][j][k][4][3] = - c2 * ( u[i][j][k][3]*u[i][j][k][1] )
	  * tmp2;
	fjac[i][j][k][4][4] = c1 * ( u[i][j][k][1] * tmp1 );
	njac[i][j][k][0][0] = 0.0;
	njac[i][j][k][0][1] = 0.0;
	njac[i][j][k][0][2] = 0.0;
	njac[i][j][k][0][3] = 0.0;
	njac[i][j][k][0][4] = 0.0;
	njac[i][j][k][1][0] = - con43 * c3c4 * tmp2 * u[i][j][k][1];
	njac[i][j][k][1][1] =   con43 * c3c4 * tmp1;
	njac[i][j][k][1][2] =   0.0;
	njac[i][j][k][1][3] =   0.0;
	njac[i][j][k][1][4] =   0.0;
	njac[i][j][k][2][0] = - c3c4 * tmp2 * u[i][j][k][2];
	njac[i][j][k][2][1] =   0.0;
	njac[i][j][k][2][2] =   c3c4 * tmp1;
	njac[i][j][k][2][3] =   0.0;
	njac[i][j][k][2][4] =   0.0;
	njac[i][j][k][3][0] = - c3c4 * tmp2 * u[i][j][k][3];
	njac[i][j][k][3][1] =   0.0;
	njac[i][j][k][3][2] =   0.0;
	njac[i][j][k][3][3] =   c3c4 * tmp1;
	njac[i][j][k][3][4] =   0.0;
	njac[i][j][k][4][0] = - ( con43 * c3c4
	  - c1345 ) * tmp3 * (pow2(u[i][j][k][1]))
	  - ( c3c4 - c1345 ) * tmp3 * (pow2(u[i][j][k][2]))
	  - ( c3c4 - c1345 ) * tmp3 * (pow2(u[i][j][k][3]))
	  - c1345 * tmp2 * u[i][j][k][4];
	njac[i][j][k][4][1] = ( con43 * c3c4
				- c1345 ) * tmp2 * u[i][j][k][1];
	njac[i][j][k][4][2] = ( c3c4 - c1345 ) * tmp2 * u[i][j][k][2];
	njac[i][j][k][4][3] = ( c3c4 - c1345 ) * tmp2 * u[i][j][k][3];
	njac[i][j][k][4][4] = ( c1345 ) * tmp1;
      }
      for (i = 1; i < grid_points[0]-1; i++) {
	tmp1 = dt * tx1;
	tmp2 = dt * tx2;
	lhs[i][j][k][AA][0][0] = - tmp2 * fjac[i-1][j][k][0][0]
	  - tmp1 * njac[i-1][j][k][0][0]
	  - tmp1 * dx1;
	lhs[i][j][k][AA][0][1] = - tmp2 * fjac[i-1][j][k][0][1]
	  - tmp1 * njac[i-1][j][k][0][1];
	lhs[i][j][k][AA][0][2] = - tmp2 * fjac[i-1][j][k][0][2]
	  - tmp1 * njac[i-1][j][k][0][2];
	lhs[i][j][k][AA][0][3] = - tmp2 * fjac[i-1][j][k][0][3]
	  - tmp1 * njac[i-1][j][k][0][3];
	lhs[i][j][k][AA][0][4] = - tmp2 * fjac[i-1][j][k][0][4]
	  - tmp1 * njac[i-1][j][k][0][4];
	lhs[i][j][k][AA][1][0] = - tmp2 * fjac[i-1][j][k][1][0]
	  - tmp1 * njac[i-1][j][k][1][0];
	lhs[i][j][k][AA][1][1] = - tmp2 * fjac[i-1][j][k][1][1]
	  - tmp1 * njac[i-1][j][k][1][1]
	  - tmp1 * dx2;
	lhs[i][j][k][AA][1][2] = - tmp2 * fjac[i-1][j][k][1][2]
	  - tmp1 * njac[i-1][j][k][1][2];
	lhs[i][j][k][AA][1][3] = - tmp2 * fjac[i-1][j][k][1][3]
	  - tmp1 * njac[i-1][j][k][1][3];
	lhs[i][j][k][AA][1][4] = - tmp2 * fjac[i-1][j][k][1][4]
	  - tmp1 * njac[i-1][j][k][1][4];
	lhs[i][j][k][AA][2][0] = - tmp2 * fjac[i-1][j][k][2][0]
	  - tmp1 * njac[i-1][j][k][2][0];
	lhs[i][j][k][AA][2][1] = - tmp2 * fjac[i-1][j][k][2][1]
	  - tmp1 * njac[i-1][j][k][2][1];
	lhs[i][j][k][AA][2][2] = - tmp2 * fjac[i-1][j][k][2][2]
	  - tmp1 * njac[i-1][j][k][2][2]
	  - tmp1 * dx3;
	lhs[i][j][k][AA][2][3] = - tmp2 * fjac[i-1][j][k][2][3]
	  - tmp1 * njac[i-1][j][k][2][3];
	lhs[i][j][k][AA][2][4] = - tmp2 * fjac[i-1][j][k][2][4]
	  - tmp1 * njac[i-1][j][k][2][4];
	lhs[i][j][k][AA][3][0] = - tmp2 * fjac[i-1][j][k][3][0]
	  - tmp1 * njac[i-1][j][k][3][0];
	lhs[i][j][k][AA][3][1] = - tmp2 * fjac[i-1][j][k][3][1]
	  - tmp1 * njac[i-1][j][k][3][1];
	lhs[i][j][k][AA][3][2] = - tmp2 * fjac[i-1][j][k][3][2]
	  - tmp1 * njac[i-1][j][k][3][2];
	lhs[i][j][k][AA][3][3] = - tmp2 * fjac[i-1][j][k][3][3]
	  - tmp1 * njac[i-1][j][k][3][3]
	  - tmp1 * dx4;
	lhs[i][j][k][AA][3][4] = - tmp2 * fjac[i-1][j][k][3][4]
	  - tmp1 * njac[i-1][j][k][3][4];
	lhs[i][j][k][AA][4][0] = - tmp2 * fjac[i-1][j][k][4][0]
	  - tmp1 * njac[i-1][j][k][4][0];
	lhs[i][j][k][AA][4][1] = - tmp2 * fjac[i-1][j][k][4][1]
	  - tmp1 * njac[i-1][j][k][4][1];
	lhs[i][j][k][AA][4][2] = - tmp2 * fjac[i-1][j][k][4][2]
	  - tmp1 * njac[i-1][j][k][4][2];
	lhs[i][j][k][AA][4][3] = - tmp2 * fjac[i-1][j][k][4][3]
	  - tmp1 * njac[i-1][j][k][4][3];
	lhs[i][j][k][AA][4][4] = - tmp2 * fjac[i-1][j][k][4][4]
	  - tmp1 * njac[i-1][j][k][4][4]
	  - tmp1 * dx5;
	lhs[i][j][k][BB][0][0] = 1.0
	  + tmp1 * 2.0 * njac[i][j][k][0][0]
	  + tmp1 * 2.0 * dx1;
	lhs[i][j][k][BB][0][1] = tmp1 * 2.0 * njac[i][j][k][0][1];
	lhs[i][j][k][BB][0][2] = tmp1 * 2.0 * njac[i][j][k][0][2];
	lhs[i][j][k][BB][0][3] = tmp1 * 2.0 * njac[i][j][k][0][3];
	lhs[i][j][k][BB][0][4] = tmp1 * 2.0 * njac[i][j][k][0][4];
	lhs[i][j][k][BB][1][0] = tmp1 * 2.0 * njac[i][j][k][1][0];
	lhs[i][j][k][BB][1][1] = 1.0
	  + tmp1 * 2.0 * njac[i][j][k][1][1]
	  + tmp1 * 2.0 * dx2;
	lhs[i][j][k][BB][1][2] = tmp1 * 2.0 * njac[i][j][k][1][2];
	lhs[i][j][k][BB][1][3] = tmp1 * 2.0 * njac[i][j][k][1][3];
	lhs[i][j][k][BB][1][4] = tmp1 * 2.0 * njac[i][j][k][1][4];
	lhs[i][j][k][BB][2][0] = tmp1 * 2.0 * njac[i][j][k][2][0];
	lhs[i][j][k][BB][2][1] = tmp1 * 2.0 * njac[i][j][k][2][1];
	lhs[i][j][k][BB][2][2] = 1.0
	  + tmp1 * 2.0 * njac[i][j][k][2][2]
	  + tmp1 * 2.0 * dx3;
	lhs[i][j][k][BB][2][3] = tmp1 * 2.0 * njac[i][j][k][2][3];
	lhs[i][j][k][BB][2][4] = tmp1 * 2.0 * njac[i][j][k][2][4];
	lhs[i][j][k][BB][3][0] = tmp1 * 2.0 * njac[i][j][k][3][0];
	lhs[i][j][k][BB][3][1] = tmp1 * 2.0 * njac[i][j][k][3][1];
	lhs[i][j][k][BB][3][2] = tmp1 * 2.0 * njac[i][j][k][3][2];
	lhs[i][j][k][BB][3][3] = 1.0
	  + tmp1 * 2.0 * njac[i][j][k][3][3]
	  + tmp1 * 2.0 * dx4;
	lhs[i][j][k][BB][3][4] = tmp1 * 2.0 * njac[i][j][k][3][4];
	lhs[i][j][k][BB][4][0] = tmp1 * 2.0 * njac[i][j][k][4][0];
	lhs[i][j][k][BB][4][1] = tmp1 * 2.0 * njac[i][j][k][4][1];
	lhs[i][j][k][BB][4][2] = tmp1 * 2.0 * njac[i][j][k][4][2];
	lhs[i][j][k][BB][4][3] = tmp1 * 2.0 * njac[i][j][k][4][3];
	lhs[i][j][k][BB][4][4] = 1.0
	  + tmp1 * 2.0 * njac[i][j][k][4][4]
	  + tmp1 * 2.0 * dx5;
	lhs[i][j][k][CC][0][0] =  tmp2 * fjac[i+1][j][k][0][0]
	  - tmp1 * njac[i+1][j][k][0][0]
	  - tmp1 * dx1;
	lhs[i][j][k][CC][0][1] =  tmp2 * fjac[i+1][j][k][0][1]
	  - tmp1 * njac[i+1][j][k][0][1];
	lhs[i][j][k][CC][0][2] =  tmp2 * fjac[i+1][j][k][0][2]
	  - tmp1 * njac[i+1][j][k][0][2];
	lhs[i][j][k][CC][0][3] =  tmp2 * fjac[i+1][j][k][0][3]
	  - tmp1 * njac[i+1][j][k][0][3];
	lhs[i][j][k][CC][0][4] =  tmp2 * fjac[i+1][j][k][0][4]
	  - tmp1 * njac[i+1][j][k][0][4];
	lhs[i][j][k][CC][1][0] =  tmp2 * fjac[i+1][j][k][1][0]
	  - tmp1 * njac[i+1][j][k][1][0];
	lhs[i][j][k][CC][1][1] =  tmp2 * fjac[i+1][j][k][1][1]
	  - tmp1 * njac[i+1][j][k][1][1]
	  - tmp1 * dx2;
	lhs[i][j][k][CC][1][2] =  tmp2 * fjac[i+1][j][k][1][2]
	  - tmp1 * njac[i+1][j][k][1][2];
	lhs[i][j][k][CC][1][3] =  tmp2 * fjac[i+1][j][k][1][3]
	  - tmp1 * njac[i+1][j][k][1][3];
	lhs[i][j][k][CC][1][4] =  tmp2 * fjac[i+1][j][k][1][4]
	  - tmp1 * njac[i+1][j][k][1][4];
	lhs[i][j][k][CC][2][0] =  tmp2 * fjac[i+1][j][k][2][0]
	  - tmp1 * njac[i+1][j][k][2][0];
	lhs[i][j][k][CC][2][1] =  tmp2 * fjac[i+1][j][k][2][1]
	  - tmp1 * njac[i+1][j][k][2][1];
	lhs[i][j][k][CC][2][2] =  tmp2 * fjac[i+1][j][k][2][2]
	  - tmp1 * njac[i+1][j][k][2][2]
	  - tmp1 * dx3;
	lhs[i][j][k][CC][2][3] =  tmp2 * fjac[i+1][j][k][2][3]
	  - tmp1 * njac[i+1][j][k][2][3];
	lhs[i][j][k][CC][2][4] =  tmp2 * fjac[i+1][j][k][2][4]
	  - tmp1 * njac[i+1][j][k][2][4];
	lhs[i][j][k][CC][3][0] =  tmp2 * fjac[i+1][j][k][3][0]
	  - tmp1 * njac[i+1][j][k][3][0];
	lhs[i][j][k][CC][3][1] =  tmp2 * fjac[i+1][j][k][3][1]
	  - tmp1 * njac[i+1][j][k][3][1];
	lhs[i][j][k][CC][3][2] =  tmp2 * fjac[i+1][j][k][3][2]
	  - tmp1 * njac[i+1][j][k][3][2];
	lhs[i][j][k][CC][3][3] =  tmp2 * fjac[i+1][j][k][3][3]
	  - tmp1 * njac[i+1][j][k][3][3]
	  - tmp1 * dx4;
	lhs[i][j][k][CC][3][4] =  tmp2 * fjac[i+1][j][k][3][4]
	  - tmp1 * njac[i+1][j][k][3][4];
	lhs[i][j][k][CC][4][0] =  tmp2 * fjac[i+1][j][k][4][0]
	  - tmp1 * njac[i+1][j][k][4][0];
	lhs[i][j][k][CC][4][1] =  tmp2 * fjac[i+1][j][k][4][1]
	  - tmp1 * njac[i+1][j][k][4][1];
	lhs[i][j][k][CC][4][2] =  tmp2 * fjac[i+1][j][k][4][2]
	  - tmp1 * njac[i+1][j][k][4][2];
	lhs[i][j][k][CC][4][3] =  tmp2 * fjac[i+1][j][k][4][3]
	  - tmp1 * njac[i+1][j][k][4][3];
	lhs[i][j][k][CC][4][4] =  tmp2 * fjac[i+1][j][k][4][4]
	  - tmp1 * njac[i+1][j][k][4][4]
	  - tmp1 * dx5;
      }
    }
  }
}