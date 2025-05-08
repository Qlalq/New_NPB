#include <stdio.h> // Assuming stdio.h is needed for completeness, though not used in the snippet
#include <math.h> // Assuming pow2 might be related to math functions, replaced with x*x

// Assuming constants like C1, C2, C3, C4, C5, dt, tx1, tx2, ty1, ty2, tz1, tz2,
// dx1..5, dy1..5, dz1..5, ist, iend, jst, jend are defined elsewhere in scope.
// Assuming arrays u, d, a, b, c are defined and accessible with appropriate dimensions.
// The dimensions are implicitly u[?][?][?][5], d[?][?][?][5][5], a[?][?][?][5][5], b[?][?][?][5][5], c[?][?][?][5][5].
// The indices used suggest d, a, b, c are indexed as [i][j][p][q] where p,q are 0-4,
// and u is indexed as [i_offset][j_offset][k_offset][var] where var is 0-4.

// Helper macro/function for square, replacing pow2
#define square(x) ((x)*(x))

static void jacu(int k) {

/*--------------------------------------------------------------------
c   compute the upper triangular part of the jacobian matrix for a given k-plane
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables used within each (i,j) iteration
--------------------------------------------------------------------*/
  int i, j;
  double  r43;
  double  c1345;
  double  c34;
  // tmp1, tmp2, tmp3 are temporary variables recalculated multiple times
  // within each (i,j) iteration, based on different u array accesses.
  // They are local to the (i,j) iteration block and do not create dependencies.
  double  tmp1, tmp2, tmp3;

  r43 = ( 4.0 / 3.0 );
  c1345 = C1 * C3 * C4 * C5;
  c34 = C3 * C4;

  // The following nested loops iterate over a 2D domain (i, j).
  // Each iteration (i, j) calculates a 5x5 block for matrices d, a, b, and c
  // corresponding to the grid point (i, j) in the k-plane.
  // The calculations for each (i, j) block are independent of other (i', j') blocks.
  // Data Dependencies:
  // - Writes to d[i][j], a[i][j], b[i][j], c[i][j] are unique for each (i,j) pair.
  // - Reads from u use indices like u[i][j][k], u[i+1][j][k], u[i][j+1][k], u[i][j][k+1].
  //   These reads are consistent for a given (i,j) block calculation and do not
  //   create loop-carried dependencies across i or j iterations.
  // Load Imbalance:
  // - The amount of work inside the inner loop is roughly constant for all (i,j).
  //   Parallelizing either the outer loop (i), the inner loop (j), or collapsing
  //   both loops would distribute work evenly assuming a sufficiently large domain.
  // Synchronization Overhead:
  // - Since writes are to distinct memory locations for each (i,j), no synchronization
  //   (like locks or atomics) is needed for accessing d, a, b, c.
  // - u is read-only.
  // - Local variables i, j, tmp1, tmp2, tmp3 are naturally private to each thread/iteration.
  // This structure is highly suitable for OpenMP parallelization on the i and/or j loops.

#if defined(_OPENMP)
  // Original OpenMP block with reversed loops. OpenMP pragmas handle parallel execution.
  #pragma omp for nowait schedule(static) private(i, j, tmp1, tmp2, tmp3) // Explicitly list private variables
  for (i = iend; i >= ist; i--) {
      for (j = jend; j >= jst; j--) {
#else
  // Standard sequential loops. This structure is parallelizable.
  // OpenMP parallelization would typically be applied to the outer loop (i)
  // or both loops (collapse).
  for (i = ist; i <= iend; i++) {
    for (j = jst; j <= jend; j++) {
#endif

/*--------------------------------------------------------------------
c   form the block daigonal d[i][j][p][q]
--------------------------------------------------------------------*/
      // Calculate temporary values based on u[i][j][k][0] for d[i][j] block
      tmp1 = 1.0 / u[i][j][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      d[i][j][0][0] =  1.0
	+ dt * 2.0 * (   tx1 * dx1
			 + ty1 * dy1
			 + tz1 * dz1 );
      d[i][j][0][1] =  0.0;
      d[i][j][0][2] =  0.0;
      d[i][j][0][3] =  0.0;
      d[i][j][0][4] =  0.0;

      d[i][j][1][0] =  dt * 2.0
	* (  tx1 * ( - r43 * c34 * tmp2 * u[i][j][k][1] )
	     + ty1 * ( -       c34 * tmp2 * u[i][j][k][1] )
	     + tz1 * ( -       c34 * tmp2 * u[i][j][k][1] ) );
      d[i][j][1][1] =  1.0
	+ dt * 2.0
	* (  tx1 * r43 * c34 * tmp1
	     + ty1 *       c34 * tmp1
	     + tz1 *       c34 * tmp1 )
	+ dt * 2.0 * (   tx1 * dx2
			 + ty1 * dy2
			 + tz1 * dz2  );
      d[i][j][1][2] = 0.0;
      d[i][j][1][3] = 0.0;
      d[i][j][1][4] = 0.0;

      d[i][j][2][0] = dt * 2.0
	* (  tx1 * ( -       c34 * tmp2 * u[i][j][k][2] )
	     + ty1 * ( - r43 * c34 * tmp2 * u[i][j][k][2] )
	     + tz1 * ( -       c34 * tmp2 * u[i][j][k][2] ) );
      d[i][j][2][1] = 0.0;
      d[i][j][2][2] = 1.0
	+ dt * 2.0
	* (  tx1 *       c34 * tmp1
	     + ty1 * r43 * c34 * tmp1
	     + tz1 *       c34 * tmp1 )
	+ dt * 2.0 * (  tx1 * dx3
			+ ty1 * dy3
			+ tz1 * dz3 );
      d[i][j][2][3] = 0.0;
      d[i][j][2][4] = 0.0;

      d[i][j][3][0] = dt * 2.0
	* (  tx1 * ( -       c34 * tmp2 * u[i][j][k][3] )
	     + ty1 * ( -       c34 * tmp2 * u[i][j][k][3] )
	     + tz1 * ( - r43 * c34 * tmp2 * u[i][j][k][3] ) );
      d[i][j][3][1] = 0.0;
      d[i][j][3][2] = 0.0;
      d[i][j][3][3] = 1.0
	+ dt * 2.0
	* (  tx1 *       c34 * tmp1
	     + ty1 *       c34 * tmp1
	     + tz1 * r43 * c34 * tmp1 )
	+ dt * 2.0 * (  tx1 * dx4
			+ ty1 * dy4
			+ tz1 * dz4 );
      d[i][j][3][4] = 0.0;

      d[i][j][4][0] = dt * 2.0
	* ( tx1 * ( - ( r43*c34 - c1345 ) * tmp3 * ( square(u[i][j][k][1]) )
		    - ( c34 - c1345 ) * tmp3 * ( square(u[i][j][k][2]) )
		    - ( c34 - c1345 ) * tmp3 * ( square(u[i][j][k][3]) )
		    - ( c1345 ) * tmp2 * u[i][j][k][4] )
	    + ty1 * ( - ( c34 - c1345 ) * tmp3 * ( square(u[i][j][k][1]) )
		      - ( r43*c34 - c1345 ) * tmp3 * ( square(u[i][j][k][2]) )
		      - ( c34 - c1345 ) * tmp3 * ( square(u[i][j][k][3]) )
		      - ( c1345 ) * tmp2 * u[i][j][k][4] )
	    + tz1 * ( - ( c34 - c1345 ) * tmp3 * ( square(u[i][j][k][1]) )
		      - ( c34 - c1345 ) * tmp3 * ( square(u[i][j][k][2]) )
		      - ( r43*c34 - c1345 ) * tmp3 * ( square(u[i][j][k][3]) )
		      - ( c1345 ) * tmp2 * u[i][j][k][4] ) );
      d[i][j][4][1] = dt * 2.0
	* ( tx1 * ( r43*c34 - c1345 ) * tmp2 * u[i][j][k][1]
	    + ty1 * (     c34 - c1345 ) * tmp2 * u[i][j][k][1]
	    + tz1 * (     c34 - c1345 ) * tmp2 * u[i][j][k][1] );
      d[i][j][4][2] = dt * 2.0
	* ( tx1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][2]
	    + ty1 * ( r43*c34 -c1345 ) * tmp2 * u[i][j][k][2]
	    + tz1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][2] );
      d[i][j][4][3] = dt * 2.0
	* ( tx1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][3]
	    + ty1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][3]
	    + tz1 * ( r43*c34 - c1345 ) * tmp2 * u[i][j][k][3] );
      d[i][j][4][4] = 1.0
        + dt * 2.0 * ( tx1 * c1345 * tmp1
		       + ty1 * c1345 * tmp1
		       + tz1 * c1345 * tmp1 )
        + dt * 2.0 * (  tx1 * dx5
			+  ty1 * dy5
			+  tz1 * dz5 );

/*--------------------------------------------------------------------
c   form the first block sub-diagonal a[i][j][p][q]
--------------------------------------------------------------------*/
      // Calculate temporary values based on u[i+1][j][k][0] for a[i][j] block
      tmp1 = 1.0 / u[i+1][j][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      a[i][j][0][0] = - dt * tx1 * dx1;
      a[i][j][0][1] =   dt * tx2;
      a[i][j][0][2] =   0.0;
      a[i][j][0][3] =   0.0;
      a[i][j][0][4] =   0.0;

      a[i][j][1][0] =  dt * tx2
	* ( - ( u[i+1][j][k][1] * tmp1 ) *( u[i+1][j][k][1] * tmp1 )
	    + C2 * 0.50 * (  u[i+1][j][k][1] * u[i+1][j][k][1]
                             + u[i+1][j][k][2] * u[i+1][j][k][2]
                             + u[i+1][j][k][3] * u[i+1][j][k][3] ) * tmp2 )
	- dt * tx1 * ( - r43 * c34 * tmp2 * u[i+1][j][k][1] );
      a[i][j][1][1] =  dt * tx2
	* ( ( 2.0 - C2 ) * ( u[i+1][j][k][1] * tmp1 ) )
	- dt * tx1 * ( r43 * c34 * tmp1 )
	- dt * tx1 * dx2;
      a[i][j][1][2] =  dt * tx2
	* ( - C2 * ( u[i+1][j][k][2] * tmp1 ) );
      a[i][j][1][3] =  dt * tx2
	* ( - C2 * ( u[i+1][j][k][3] * tmp1 ) );
      a[i][j][1][4] =  dt * tx2 * C2 ;

      a[i][j][2][0] =  dt * tx2
	* ( - ( u[i+1][j][k][1] * u[i+1][j][k][2] ) * tmp2 )
	- dt * tx1 * ( - c34 * tmp2 * u[i+1][j][k][2] );
      a[i][j][2][1] =  dt * tx2 * ( u[i+1][j][k][2] * tmp1 );
      a[i][j][2][2] =  dt * tx2 * ( u[i+1][j][k][1] * tmp1 )
	- dt * tx1 * ( c34 * tmp1 )
	- dt * tx1 * dx3;
      a[i][j][2][3] = 0.0;
      a[i][j][2][4] = 0.0;

      a[i][j][3][0] = dt * tx2
	* ( - ( u[i+1][j][k][1]*u[i+1][j][k][3] ) * tmp2 )
	- dt * tx1 * ( - c34 * tmp2 * u[i+1][j][k][3] );
      a[i][j][3][1] = dt * tx2 * ( u[i+1][j][k][3] * tmp1 );
      a[i][j][3][2] = 0.0;
      a[i][j][3][3] = dt * tx2 * ( u[i+1][j][k][1] * tmp1 )
	- dt * tx1 * ( c34 * tmp1 )
	- dt * tx1 * dx4;
      a[i][j][3][4] = 0.0;

      a[i][j][4][0] = dt * tx2
	* ( ( C2 * (  u[i+1][j][k][1] * u[i+1][j][k][1]
		      + u[i+1][j][k][2] * u[i+1][j][k][2]
		      + u[i+1][j][k][3] * u[i+1][j][k][3] ) * tmp2
	      - C1 * ( u[i+1][j][k][4] * tmp1 ) )
	    * ( u[i+1][j][k][1] * tmp1 ) )
	- dt * tx1
	* ( - ( r43*c34 - c1345 ) * tmp3 * ( square(u[i+1][j][k][1]) )
	    - (     c34 - c1345 ) * tmp3 * ( square(u[i+1][j][k][2]) )
	    - (     c34 - c1345 ) * tmp3 * ( square(u[i+1][j][k][3]) )
	    - c1345 * tmp2 * u[i+1][j][k][4] );
      a[i][j][4][1] = dt * tx2
	* ( C1 * ( u[i+1][j][k][4] * tmp1 )
	    - 0.50 * C2
	    * ( (  3.0*u[i+1][j][k][1]*u[i+1][j][k][1]
		   + u[i+1][j][k][2]*u[i+1][j][k][2]
		   + u[i+1][j][k][3]*u[i+1][j][k][3] ) * tmp2 ) )
	- dt * tx1
	* ( r43*c34 - c1345 ) * tmp2 * u[i+1][j][k][1];
      a[i][j][4][2] = dt * tx2
	* ( - C2 * ( u[i+1][j][k][2]*u[i+1][j][k][1] ) * tmp2 )
	- dt * tx1
	* (  c34 - c1345 ) * tmp2 * u[i+1][j][k][2];
      a[i][j][4][3] = dt * tx2
	* ( - C2 * ( u[i+1][j][k][3]*u[i+1][j][k][1] ) * tmp2 )
	- dt * tx1
	* (  c34 - c1345 ) * tmp2 * u[i+1][j][k][3];
      a[i][j][4][4] = dt * tx2
	* ( C1 * ( u[i+1][j][k][1] * tmp1 ) )
	- dt * tx1 * c1345 * tmp1
	- dt * tx1 * dx5;

/*--------------------------------------------------------------------
c   form the second block sub-diagonal b[i][j][p][q]
--------------------------------------------------------------------*/
      // Calculate temporary values based on u[i][j+1][k][0] for b[i][j] block
      tmp1 = 1.0 / u[i][j+1][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      b[i][j][0][0] = - dt * ty1 * dy1;
      b[i][j][0][1] =   0.0;
      b[i][j][0][2] =  dt * ty2;
      b[i][j][0][3] =   0.0;
      b[i][j][0][4] =   0.0;

      b[i][j][1][0] =  dt * ty2
	* ( - ( u[i][j+1][k][1]*u[i][j+1][k][2] ) * tmp2 )
	- dt * ty1 * ( - c34 * tmp2 * u[i][j+1][k][1] );
      b[i][j][1][1] =  dt * ty2 * ( u[i][j+1][k][2] * tmp1 )
	- dt * ty1 * ( c34 * tmp1 )
	- dt * ty1 * dy2;
      b[i][j][1][2] =  dt * ty2 * ( u[i][j+1][k][1] * tmp1 );
      b[i][j][1][3] = 0.0;
      b[i][j][1][4] = 0.0;

      b[i][j][2][0] =  dt * ty2
	* ( - ( u[i][j+1][k][2] * tmp1 ) *( u[i][j+1][k][2] * tmp1 )
	    + 0.50 * C2 * ( (  u[i][j+1][k][1] * u[i][j+1][k][1]
			       + u[i][j+1][k][2] * u[i][j+1][k][2]
			       + u[i][j+1][k][3] * u[i][j+1][k][3] )
			    * tmp2 ) )
	- dt * ty1 * ( - r43 * c34 * tmp2 * u[i][j+1][k][2] );
      b[i][j][2][1] =  dt * ty2
	* ( - C2 * ( u[i][j+1][k][1] * tmp1 ) );
      b[i][j][2][2] =  dt * ty2 * ( ( 2.0 - C2 )
				 * ( u[i][j+1][k][2] * tmp1 ) )
	- dt * ty1 * ( r43 * c34 * tmp1 )
	- dt * ty1 * dy3;
      b[i][j][2][3] =  dt * ty2
	* ( - C2 * ( u[i][j+1][k][3] * tmp1 ) );
      b[i][j][2][4] =  dt * ty2 * C2;

      b[i][j][3][0] =  dt * ty2
	* ( - ( u[i][j+1][k][2]*u[i][j+1][k][3] ) * tmp2 )
	- dt * ty1 * ( - c34 * tmp2 * u[i][j+1][k][3] );
      b[i][j][3][1] = 0.0;
      b[i][j][3][2] =  dt * ty2 * ( u[i][j+1][k][3] * tmp1 );
      b[i][j][3][3] =  dt * ty2 * ( u[i][j+1][k][2] * tmp1 )
	- dt * ty1 * ( c34 * tmp1 )
	- dt * ty1 * dy4;
      b[i][j][3][4] = 0.0;

      b[i][j][4][0] =  dt * ty2
	* ( ( C2 * (  u[i][j+1][k][1] * u[i][j+1][k][1]
		      + u[i][j+1][k][2] * u[i][j+1][k][2]
		      + u[i][j+1][k][3] * u[i][j+1][k][3] ) * tmp2
	      - C1 * ( u[i][j+1][k][4] * tmp1 ) )
	    * ( u[i][j+1][k][2] * tmp1 ) )
	- dt * ty1
	* ( - (     c34 - c1345 )*tmp3*( square(u[i][j+1][k][1]) )
	    - ( r43*c34 - c1345 )*tmp3*( square(u[i][j+1][k][2]) )
	    - (     c34 - c1345 )*tmp3*( square(u[i][j+1][k][3]) )
	    - c1345*tmp2*u[i][j+1][k][4] );
      b[i][j][4][1] =  dt * ty2
	* ( - C2 * ( u[i][j+1][k][1]*u[i][j+1][k][2] ) * tmp2 )
	- dt * ty1
	* ( c34 - c1345 ) * tmp2 * u[i][j+1][k][1];
      b[i][j][4][2] =  dt * ty2
	* ( C1 * ( u[i][j+1][k][4] * tmp1 )
	    - 0.50 * C2
	    * ( (  u[i][j+1][k][1]*u[i][j+1][k][1]
		   + 3.0 * u[i][j+1][k][2]*u[i][j+1][k][2]
		   + u[i][j+1][k][3]*u[i][j+1][k][3] ) * tmp2 ) )
	- dt * ty1
	* ( r43*c34 - c1345 ) * tmp2 * u[i][j+1][k][2];
      b[i][j][4][3] =  dt * ty2
	* ( - C2 * ( u[i][j+1][k][2]*u[i][j+1][k][3] ) * tmp2 )
	- dt * ty1 * ( c34 - c1345 ) * tmp2 * u[i][j+1][k][3];
      b[i][j][4][4] =  dt * ty2
	* ( C1 * ( u[i][j+1][k][2] * tmp1 ) )
	- dt * ty1 * c1345 * tmp1
	- dt * ty1 * dy5;

/*--------------------------------------------------------------------
c   form the third block sub-diagonal c[i][j][p][q]
--------------------------------------------------------------------*/
      // Calculate temporary values based on u[i][j][k+1][0] for c[i][j] block
      tmp1 = 1.0 / u[i][j][k+1][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      c[i][j][0][0] = - dt * tz1 * dz1;
      c[i][j][0][1] =   0.0;
      c[i][j][0][2] =   0.0;
      c[i][j][0][3] = dt * tz2;
      c[i][j][0][4] =   0.0;

      c[i][j][1][0] = dt * tz2
	* ( - ( u[i][j][k+1][1]*u[i][j][k+1][3] ) * tmp2 )
	- dt * tz1 * ( - c34 * tmp2 * u[i][j][k+1][1] );
      c[i][j][1][1] = dt * tz2 * ( u[i][j][k+1][3] * tmp1 )
	- dt * tz1 * c34 * tmp1
	- dt * tz1 * dz2 ;
      c[i][j][1][2] = 0.0;
      c[i][j][1][3] = dt * tz2 * ( u[i][j][k+1][1] * tmp1 );
      c[i][j][1][4] = 0.0;

      c[i][j][2][0] = dt * tz2
	* ( - ( u[i][j][k+1][2]*u[i][j][k+1][3] ) * tmp2 )
	- dt * tz1 * ( - c34 * tmp2 * u[i][j][k+1][2] );
      c[i][j][2][1] = 0.0;
      c[i][j][2][2] = dt * tz2 * ( u[i][j][k+1][3] * tmp1 )
	- dt * tz1 * ( c34 * tmp1 )
	- dt * tz1 * dz3;
      c[i][j][2][3] = dt * tz2 * ( u[i][j][k+1][2] * tmp1 );
      c[i][j][2][4] = 0.0;

      c[i][j][3][0] = dt * tz2
	* ( - ( u[i][j][k+1][3] * tmp1 ) *( u[i][j][k+1][3] * tmp1 )
	    + 0.50 * C2
	    * ( ( u[i][j][k+1][1] * u[i][j][k+1][1]
		  + u[i][j][k+1][2] * u[i][j][k+1][2]
		  + u[i][j][k+1][3] * u[i][j][k+1][3] ) * tmp2 ) )
	- dt * tz1 * ( - r43 * c34 * tmp2 * u[i][j][k+1][3] );
      c[i][j][3][1] = dt * tz2
	* ( - C2 * ( u[i][j][k+1][1] * tmp1 ) );
      c[i][j][3][2] = dt * tz2
	* ( - C2 * ( u[i][j][k+1][2] * tmp1 ) );
      c[i][j][3][3] = dt * tz2 * ( 2.0 - C2 )
	* ( u[i][j][k+1][3] * tmp1 )
	- dt * tz1 * ( r43 * c34 * tmp1 )
	- dt * tz1 * dz4;
      c[i][j][3][4] = dt * tz2 * C2;

      c[i][j][4][0] = dt * tz2
	* ( ( C2 * (  u[i][j][k+1][1] * u[i][j][k+1][1]
                      + u[i][j][k+1][2] * u[i][j][k+1][2]
                      + u[i][j][k+1][3] * u[i][j][k+1][3] ) * tmp2
	      - C1 * ( u[i][j][k+1][4] * tmp1 ) )
	    * ( u[i][j][k+1][3] * tmp1 ) )
	- dt * tz1
	* ( - ( c34 - c1345 ) * tmp3 * ( square(u[i][j][k+1][1]) )
	    - ( c34 - c1345 ) * tmp3 * ( square(u[i][j][k+1][2]) )
	    - ( r43*c34 - c1345 )* tmp3 * ( square(u[i][j][k+1][3]) )
	    - c1345 * tmp2 * u[i][j][k+1][4] );
      c[i][j][4][1] = dt * tz2
	* ( - C2 * ( u[i][j][k+1][1]*u[i][j][k+1][3] ) * tmp2 )
	- dt * tz1 * ( c34 - c1345 ) * tmp2 * u[i][j][k+1][1];
      c[i][j][4][2] = dt * tz2
	* ( - C2 * ( u[i][j][k+1][2]*u[i][j][k+1][3] ) * tmp2 )
	- dt * tz1 * ( c34 - c1345 ) * tmp2 * u[i][j][k+1][2];
      c[i][j][4][3] = dt * tz2
	* ( C1 * ( u[i][j][k+1][4] * tmp1 )
            - 0.50 * C2
            * ( (  u[i][j][k+1][1]*u[i][j][k+1][1]
		   + u[i][j][k+1][2]*u[i][j][k+1][2]
		   + 3.0*u[i][j][k+1][3]*u[i][j][k+1][3] ) * tmp2 ) )
	- dt * tz1 * ( r43*c34 - c1345 ) * tmp2 * u[i][j][k+1][3];
      c[i][j][4][4] = dt * tz2
	* ( C1 * ( u[i][j][k+1][3] * tmp1 ) )
	- dt * tz1 * c1345 * tmp1
	- dt * tz1 * dz5;
    }
  }
}