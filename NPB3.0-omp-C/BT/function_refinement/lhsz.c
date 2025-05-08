// Note: Assumed definitions for grid_points, u, fjac, njac, lhs,
// c1, c2, c3c4, con43, c1345, dt, tz1, tz2, dz1, dz2, dz3, dz4, dz5,
// AA, BB, CC, and pow2() exist outside this function, likely as global
// variables/constants or via included headers. No headers are added
// in the output as requested.

static void lhsz(void) {
  int i, j, k;
  // Declare temporary variables locally for better scoping and
  // to facilitate potential thread-private allocation in parallel regions.
  double tmp1, tmp2, tmp3;

  // This first block computes the fjac and njac arrays.
  // The calculation for each element fjac[i][j][k] and njac[i][j][k]
  // depends only on u[i][j][k] and constants.
  // There are no loop-carried dependencies across i, j, or k.
  // This loop nest is fully parallelizable.
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 0; k < grid_points[2]; k++) {
	tmp1 = 1.0 / u[i][j][k][0];
	tmp2 = tmp1 * tmp1;
	tmp3 = tmp1 * tmp2;

	fjac[i][j][k][0][0] = 0.0;
	fjac[i][j][k][0][1] = 0.0;
	fjac[i][j][k][0][2] = 0.0;
	fjac[i][j][k][0][3] = 1.0;
	fjac[i][j][k][0][4] = 0.0;
	fjac[i][j][k][1][0] = - ( u[i][j][k][1]*u[i][j][k][3] )
	  * tmp2;
	fjac[i][j][k][1][1] = u[i][j][k][3] * tmp1;
	fjac[i][j][k][1][2] = 0.0;
	fjac[i][j][k][1][3] = u[i][j][k][1] * tmp1;
	fjac[i][j][k][1][4] = 0.0;
	fjac[i][j][k][2][0] = - ( u[i][j][k][2]*u[i][j][k][3] )
	  * tmp2;
	fjac[i][j][k][2][1] = 0.0;
	fjac[i][j][k][2][2] = u[i][j][k][3] * tmp1;
	fjac[i][j][k][2][3] = u[i][j][k][2] * tmp1;
	fjac[i][j][k][2][4] = 0.0;
	fjac[i][j][k][3][0] = - (u[i][j][k][3]*u[i][j][k][3] * tmp2 )
	  + 0.50 * c2 * ( (  u[i][j][k][1] * u[i][j][k][1]
			     + u[i][j][k][2] * u[i][j][k][2]
			     + u[i][j][k][3] * u[i][j][k][3] ) * tmp2 );
	fjac[i][j][k][3][1] = - c2 *  u[i][j][k][1] * tmp1;
	fjac[i][j][k][3][2] = - c2 *  u[i][j][k][2] * tmp1;
	fjac[i][j][k][3][3] = ( 2.0 - c2 )
	  *  u[i][j][k][3] * tmp1;
	fjac[i][j][k][3][4] = c2;
	fjac[i][j][k][4][0] = ( c2 * (  u[i][j][k][1] * u[i][j][k][1]
					+ u[i][j][k][2] * u[i][j][k][2]
					+ u[i][j][k][3] * u[i][j][k][3] )
				* tmp2
				- c1 * ( u[i][j][k][4] * tmp1 ) )
	  * ( u[i][j][k][3] * tmp1 );
	fjac[i][j][k][4][1] = - c2 * ( u[i][j][k][1]*u[i][j][k][3] )
	  * tmp2;
	fjac[i][j][k][4][2] = - c2 * ( u[i][j][k][2]*u[i][j][k][3] )
	  * tmp2;
	fjac[i][j][k][4][3] = c1 * ( u[i][j][k][4] * tmp1 )
	  - 0.50 * c2
	  * ( (  u[i][j][k][1]*u[i][j][k][1]
		 + u[i][j][k][2]*u[i][j][k][2]
		 + 3.0*u[i][j][k][3]*u[i][j][k][3] )
	      * tmp2 );
	fjac[i][j][k][4][4] = c1 * u[i][j][k][3] * tmp1;

	njac[i][j][k][0][0] = 0.0;
	njac[i][j][k][0][1] = 0.0;
	njac[i][j][k][0][2] = 0.0;
	njac[i][j][k][0][3] = 0.0;
	njac[i][j][k][0][4] = 0.0;
	njac[i][j][k][1][0] = - c3c4 * tmp2 * u[i][j][k][1];
	njac[i][j][k][1][1] =   c3c4 * tmp1;
	njac[i][j][k][1][2] =   0.0;
	njac[i][j][k][1][3] =   0.0;
	njac[i][j][k][1][4] =   0.0;
	njac[i][j][k][2][0] = - c3c4 * tmp2 * u[i][j][k][2];
	njac[i][j][k][2][1] =   0.0;
	njac[i][j][k][2][2] =   c3c4 * tmp1;
	njac[i][j][k][2][3] =   0.0;
	njac[i][j][k][2][4] =   0.0;
	njac[i][j][k][3][0] = - con43 * c3c4 * tmp2 * u[i][j][k][3];
	njac[i][j][k][3][1] =   0.0;
	njac[i][j][k][3][2] =   0.0;
	njac[i][j][k][3][3] =   con43 * c3 * c4 * tmp1;
	njac[i][j][k][3][4] =   0.0;
	njac[i][j][k][4][0] = - (  c3c4
	  - c1345 ) * tmp3 * (pow2(u[i][j][k][1]))
	  - ( c3c4 - c1345 ) * tmp3 * (pow2(u[i][j][k][2]))
	  - ( con43 * c3c4
	      - c1345 ) * tmp3 * (pow2(u[i][j][k][3]))
	  - c1345 * tmp2 * u[i][j][k][4];
	njac[i][j][k][4][1] = (  c3c4 - c1345 ) * tmp2 * u[i][j][k][1];
	njac[i][j][k][4][2] = (  c3c4 - c1345 ) * tmp2 * u[i][j][k][2];
	njac[i][j][k][4][3] = ( con43 * c3c4
				- c1345 ) * tmp2 * u[i][j][k][3];
	njac[i][j][k][4][4] = ( c1345 )* tmp1;
      }
    }
  }

  // This second block computes the lhs array using the fjac and njac values.
  // The calculation for each element lhs[i][j][k] depends on fjac/njac
  // at k-1, k, and k+1 indices, but these fjac/njac values were fully computed
  // and stored in the previous block.
  // There are no loop-carried dependencies on 'lhs' across i, j, or k iterations.
  // This loop nest is also fully parallelizable.

  // Reuse temporary variables for calculations in this block.
  tmp1 = dt * tz1; // Note: tmp1 and tmp2 are reassigned
  tmp2 = dt * tz2;

  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	lhs[i][j][k][AA][0][0] = - tmp2 * fjac[i][j][k-1][0][0]
	  - tmp1 * njac[i][j][k-1][0][0]
	  - tmp1 * dz1;
	lhs[i][j][k][AA][0][1] = - tmp2 * fjac[i][j][k-1][0][1]
	  - tmp1 * njac[i][j][k-1][0][1];
	lhs[i][j][k][AA][0][2] = - tmp2 * fjac[i][j][k-1][0][2]
	  - tmp1 * njac[i][j][k-1][0][2];
	lhs[i][j][k][AA][0][3] = - tmp2 * fjac[i][j][k-1][0][3]
	  - tmp1 * njac[i][j][k-1][0][3];
	lhs[i][j][k][AA][0][4] = - tmp2 * fjac[i][j][k-1][0][4]
	  - tmp1 * njac[i][j][k-1][0][4];
	lhs[i][j][k][AA][1][0] = - tmp2 * fjac[i][j][k-1][1][0]
	  - tmp1 * njac[i][j][k-1][1][0];
	lhs[i][j][k][AA][1][1] = - tmp2 * fjac[i][j][k-1][1][1]
	  - tmp1 * njac[i][j][k-1][1][1]
	  - tmp1 * dz2;
	lhs[i][j][k][AA][1][2] = - tmp2 * fjac[i][j][k-1][1][2]
	  - tmp1 * njac[i][j][k-1][1][2];
	lhs[i][j][k][AA][1][3] = - tmp2 * fjac[i][j][k-1][1][3]
	  - tmp1 * njac[i][j][k-1][1][3];
	lhs[i][j][k][AA][1][4] = - tmp2 * fjac[i][j][k-1][1][4]
	  - tmp1 * njac[i][j][k-1][1][4];
	lhs[i][j][k][AA][2][0] = - tmp2 * fjac[i][j][k-1][2][0]
	  - tmp1 * njac[i][j][k-1][2][0];
	lhs[i][j][k][AA][2][1] = - tmp2 * fjac[i][j][k-1][2][1]
	  - tmp1 * njac[i][j][k-1][2][1];
	lhs[i][j][k][AA][2][2] = - tmp2 * fjac[i][j][k-1][2][2]
	  - tmp1 * njac[i][j][k-1][2][2]
	  - tmp1 * dz3;
	lhs[i][j][k][AA][2][3] = - tmp2 * fjac[i][j][k-1][2][3]
	  - tmp1 * njac[i][j][k-1][2][3];
	lhs[i][j][k][AA][2][4] = - tmp2 * fjac[i][j][k-1][2][4]
	  - tmp1 * njac[i][j][k-1][2][4];
	lhs[i][j][k][AA][3][0] = - tmp2 * fjac[i][j][k-1][3][0]
	  - tmp1 * njac[i][j][k-1][3][0];
	lhs[i][j][k][AA][3][1] = - tmp2 * fjac[i][j][k-1][3][1]
	  - tmp1 * njac[i][j][k-1][3][1];
	lhs[i][j][k][AA][3][2] = - tmp2 * fjac[i][j][k-1][3][2]
	  - tmp1 * njac[i][j][k-1][3][2];
	lhs[i][j][k][AA][3][3] = - tmp2 * fjac[i][j][k-1][3][3]
	  - tmp1 * njac[i][j][k-1][3][3]
	  - tmp1 * dz4;
	lhs[i][j][k][AA][3][4] = - tmp2 * fjac[i][j][k-1][3][4]
	  - tmp1 * njac[i][j][k-1][3][4];
	lhs[i][j][k][AA][4][0] = - tmp2 * fjac[i][j][k-1][4][0]
	  - tmp1 * njac[i][j][k-1][4][0];
	lhs[i][j][k][AA][4][1] = - tmp2 * fjac[i][j][k-1][4][1]
	  - tmp1 * njac[i][j][k-1][4][1];
	lhs[i][j][k][AA][4][2] = - tmp2 * fjac[i][j][k-1][4][2]
	  - tmp1 * njac[i][j][k-1][4][2];
	lhs[i][j][k][AA][4][3] = - tmp2 * fjac[i][j][k-1][4][3]
	  - tmp1 * njac[i][j][k-1][4][3];
	lhs[i][j][k][AA][4][4] = - tmp2 * fjac[i][j][k-1][4][4]
	  - tmp1 * njac[i][j][k-1][4][4]
	  - tmp1 * dz5;

	lhs[i][j][k][BB][0][0] = 1.0
	  + tmp1 * 2.0 * njac[i][j][k][0][0]
	  + tmp1 * 2.0 * dz1;
	lhs[i][j][k][BB][0][1] = tmp1 * 2.0 * njac[i][j][k][0][1];
	lhs[i][j][k][BB][0][2] = tmp1 * 2.0 * njac[i][j][k][0][2];
	lhs[i][j][k][BB][0][3] = tmp1 * 2.0 * njac[i][j][k][0][3];
	lhs[i][j][k][BB][0][4] = tmp1 * 2.0 * njac[i][j][k][0][4];
	lhs[i][j][k][BB][1][0] = tmp1 * 2.0 * njac[i][j][k][1][0];
	lhs[i][j][k][BB][1][1] = 1.0
	  + tmp1 * 2.0 * njac[i][j][k][1][1]
	  + tmp1 * 2.0 * dz2;
	lhs[i][j][k][BB][1][2] = tmp1 * 2.0 * njac[i][j][k][1][2];
	lhs[i][j][k][BB][1][3] = tmp1 * 2.0 * njac[i][j][k][1][3];
	lhs[i][j][k][BB][1][4] = tmp1 * 2.0 * njac[i][j][k][1][4];
	lhs[i][j][k][BB][2][0] = tmp1 * 2.0 * njac[i][j][k][2][0];
	lhs[i][j][k][BB][2][1] = tmp1 * 2.0 * njac[i][j][k][2][1];
	lhs[i][j][k][BB][2][2] = 1.0
	  + tmp1 * 2.0 * njac[i][j][k][2][2]
	  + tmp1 * 2.0 * dz3;
	lhs[i][j][k][BB][2][3] = tmp1 * 2.0 * njac[i][j][k][2][3];
	lhs[i][j][k][BB][2][4] = tmp1 * 2.0 * njac[i][j][k][2][4];
	lhs[i][j][k][BB][3][0] = tmp1 * 2.0 * njac[i][j][k][3][0];
	lhs[i][j][k][BB][3][1] = tmp1 * 2.0 * njac[i][j][k][3][1];
	lhs[i][j][k][BB][3][2] = tmp1 * 2.0 * njac[i][j][k][3][2];
	lhs[i][j][k][BB][3][3] = 1.0
	  + tmp1 * 2.0 * njac[i][j][k][3][3]
	  + tmp1 * 2.0 * dz4;
	lhs[i][j][k][BB][3][4] = tmp1 * 2.0 * njac[i][j][k][3][4];
	lhs[i][j][k][BB][4][0] = tmp1 * 2.0 * njac[i][j][k][4][0];
	lhs[i][j][k][BB][4][1] = tmp1 * 2.0 * njac[i][j][k][4][1];
	lhs[i][j][k][BB][4][2] = tmp1 * 2.0 * njac[i][j][k][4][2];
	lhs[i][j][k][BB][4][3] = tmp1 * 2.0 * njac[i][j][k][4][3];
	lhs[i][j][k][BB][4][4] = 1.0
	  + tmp1 * 2.0 * njac[i][j][k][4][4]
	  + tmp1 * 2.0 * dz5;

	lhs[i][j][k][CC][0][0] =  tmp2 * fjac[i][j][k+1][0][0]
	  - tmp1 * njac[i][j][k+1][0][0]
	  - tmp1 * dz1;
	lhs[i][j][k][CC][0][1] =  tmp2 * fjac[i][j][k+1][0][1]
	  - tmp1 * njac[i][j][k+1][0][1];
	lhs[i][j][k][CC][0][2] =  tmp2 * fjac[i][j][k+1][0][2]
	  - tmp1 * njac[i][j][k+1][0][2];
	lhs[i][j][k][CC][0][3] =  tmp2 * fjac[i][j][k+1][0][3]
	  - tmp1 * njac[i][j][k+1][0][3];
	lhs[i][j][k][CC][0][4] =  tmp2 * fjac[i][j][k+1][0][4]
	  - tmp1 * njac[i][j][k+1][0][4];
	lhs[i][j][k][CC][1][0] =  tmp2 * fjac[i][j][k+1][1][0]
	  - tmp1 * njac[i][j][k+1][1][0];
	lhs[i][j][k][CC][1][1] =  tmp2 * fjac[i][j][k+1][1][1]
	  - tmp1 * njac[i][j][k+1][1][1]
	  - tmp1 * dz2;
	lhs[i][j][k][CC][1][2] =  tmp2 * fjac[i][j][k+1][1][2]
	  - tmp1 * njac[i][j][k+1][1][2];
	lhs[i][j][k][CC][1][3] =  tmp2 * fjac[i][j][k+1][1][3]
	  - tmp1 * njac[i][j][k+1][1][3];
	lhs[i][j][k][CC][1][4] =  tmp2 * fjac[i][j][k+1][1][4]
	  - tmp1 * njac[i][j][k+1][1][4];
	lhs[i][j][k][CC][2][0] =  tmp2 * fjac[i][j][k+1][2][0]
	  - tmp1 * njac[i][j][k+1][2][0];
	lhs[i][j][k][CC][2][1] =  tmp2 * fjac[i][j][k+1][2][1]
	  - tmp1 * njac[i][j][k+1][2][1];
	lhs[i][j][k][CC][2][2] =  tmp2 * fjac[i][j][k+1][2][2]
	  - tmp1 * njac[i][j][k+1][2][2]
	  - tmp1 * dz3;
	lhs[i][j][k][CC][2][3] =  tmp2 * fjac[i][j][k+1][2][3]
	  - tmp1 * njac[i][j][k+1][2][3];
	lhs[i][j][k][CC][2][4] =  tmp2 * fjac[i][j][k+1][2][4]
	  - tmp1 * njac[i][j][k+1][2][4];
	lhs[i][j][k][CC][3][0] =  tmp2 * fjac[i][j][k+1][3][0]
	  - tmp1 * njac[i][j][k+1][3][0];
	lhs[i][j][k][CC][3][1] =  tmp2 * fjac[i][j][k+1][3][1]
	  - tmp1 * njac[i][j][k+1][3][1];
	lhs[i][j][k][CC][3][2] =  tmp2 * fjac[i][j][k+1][3][2]
	  - tmp1 * njac[i][j][k+1][3][2];
	lhs[i][j][k][CC][3][3] =  tmp2 * fjac[i][j][k+1][3][3]
	  - tmp1 * njac[i][j][k+1][3][3]
	  - tmp1 * dz4;
	lhs[i][j][k][CC][3][4] =  tmp2 * fjac[i][j][k+1][3][4]
	  - tmp1 * njac[i][j][k+1][3][4];
	lhs[i][j][k][CC][4][0] =  tmp2 * fjac[i][j][k+1][4][0]
	  - tmp1 * njac[i][j][k+1][4][0];
	lhs[i][j][k][CC][4][1] =  tmp2 * fjac[i][j][k+1][4][1]
	  - tmp1 * njac[i][j][k+1][4][1];
	lhs[i][j][k][CC][4][2] =  tmp2 * fjac[i][j][k+1][4][2]
	  - tmp1 * njac[i][j][k+1][4][2];
	lhs[i][j][k][CC][4][3] =  tmp2 * fjac[i][j][k+1][4][3]
	  - tmp1 * njac[i][j][k+1][4][3];
	lhs[i][j][k][CC][4][4] =  tmp2 * fjac[i][j][k+1][4][4]
	  - tmp1 * njac[i][j][k+1][4][4]
	  - tmp1 * dz5;
      }
    }
  }
}