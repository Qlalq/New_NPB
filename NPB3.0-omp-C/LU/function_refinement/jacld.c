#include <stdio.h> // Assuming stdio.h is needed for context, though not strictly in the function itself

// Assume these constants and arrays are declared and accessible globally or in an outer scope
// const double C1, C2, C3, C4, C5, dt, tx1, tx2, ty1, ty2, tz1, tz2, dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, dy5, dz1, dz2, dz3, dz4, dz5;
// int ist, iend, jst, jend;
// double u[...][...][...][...], a[...][...][...][...], b[...][...][...][...], c[...][...][...][...], d[...][...][...][...];

// Helper macro for squaring, used in original code
#ifndef pow2
#define pow2(x) ((x)*(x))
#endif

/*
 * Refined jacld function for improved structure towards efficient OpenMP parallelization.
 *
 * Refinements applied:
 * 1. Data Dependencies: Analyzed loops (i and j) and confirmed there are no loop-carried
 *    dependencies on the output arrays (a, b, c, d) that would prevent parallelization
 *    of the i and j loops. Reads from the input array (u) involve neighboring elements,
 *    but u is read-only within this function, which is safe for parallel reads.
 *    No major structural changes needed for decoupling dependencies.
 * 2. Load Imbalance: The work done inside the inner loop is computationally uniform
 *    across different (i, j) pairs. The current nested loop structure provides
 *    independent work units ((i, j) iterations) which are suitable for balanced
 *    distribution across threads using standard OpenMP loop pragmas. No structural
 *    changes needed for load balancing beyond the inherent loop structure.
 * 3. Synchronization Overhead: Writes to output arrays a, b, c, d are local to
 *    the current (i, j) iteration block. There are no writes to shared data that
 *    require explicit synchronization (like locks or atomics) if (i, j) iterations
 *    are distributed among threads. Temporary variables (tmp1, tmp2, tmp3) are
 *    specific to each (i, j) block's computation and should be thread-private
 *    when parallelized. Variables computed once before the loops are loop-invariant
 *    and can be shared or first-private.
 *
 * Structure changes implemented:
 * - Hoisted constant calculations (involving dt, tx/ty/tz, dx/dy/dz, C constants)
 *   outside the main loops to avoid redundant computations. This improves serial
 *   efficiency and doesn't negatively impact parallelization.
 * - Explicitly declared temporary variables tmp1, tmp2, tmp3 just before the inner
 *   loop to highlight their scope as private to each (i, j) iteration block when
 *   parallelized. (Note: In C89, declarations must be at the start of a block.
 *   In C99+, they can be placed closer to first use). Placing before the inner loop
 *   is a common practice suitable for OpenMP `private` clause.
 *
 * This structure is well-suited for applying OpenMP pragmas to the i loop, the j loop,
 * or both (using collapse).
 */
static void jacld(int k) {
  int i, j;

  // --- Hoisted Constants / Pre-calculations ---
  // These are loop-invariant and can be calculated once before the loops.
  double r43 = ( 4.0 / 3.0 );
  double c1345 = C1 * C3 * C4 * C5;
  double c34 = C3 * C4;

  double dt_2 = dt * 2.0;

  double dt_2_tx1 = dt_2 * tx1;
  double dt_2_ty1 = dt_2 * ty1;
  double dt_2_tz1 = dt_2 * tz1;

  double dt_tx1 = dt * tx1;
  double dt_tx2 = dt * tx2;
  double dt_ty1 = dt * ty1;
  double dt_ty2 = dt * ty2;
  double dt_tz1 = dt * tz1;
  double dt_tz2 = dt * tz2;

  // Terms involving dx, dy, dz
  double dt_2_tx1_dx1 = dt_2_tx1 * dx1;
  double dt_2_ty1_dy1 = dt_2_ty1 * dy1;
  double dt_2_tz1_dz1 = dt_2_tz1 * dz1;

  double dt_2_tx1_dx2 = dt_2_tx1 * dx2;
  double dt_2_ty1_dy2 = dt_2_ty1 * dy2;
  double dt_2_tz1_dz2 = dt_2_tz1 * dz2;

  double dt_2_tx1_dx3 = dt_2_tx1 * dx3;
  double dt_2_ty1_dy3 = dt_2_ty1 * dy3;
  double dt_2_tz1_dz3 = dt_2_tz1 * dz3;

  double dt_2_tx1_dx4 = dt_2_tx1 * dx4;
  double dt_2_ty1_dy4 = dt_2_ty1 * dy4;
  double dt_2_tz1_dz4 = dt_2_tz1 * dz4;

  double dt_2_tx1_dx5 = dt_2_tx1 * dx5;
  double dt_2_ty1_dy5 = dt_2_ty1 * dy5;
  double dt_2_tz1_dz5 = dt_2_tz1 * dz5;

  // Combinations of C constants with r43/c34
  double r43_c34 = r43 * c34; // This factor appears often with tmp1
  double c34_m_c1345 = c34 - c1345;
  double r43_c34_m_c1345 = r43_c34 - c1345;
  double c2_05 = 0.50 * C2; // Factor 0.5 * C2 appears often

  // --- Main Loops ---
  // These loops define the independent computational units (i, j)
  for (i = ist; i <= iend; i++) {
    // Declare temporary variables here. These will be private to each thread's i iteration range.
    // When parallelizing the j loop or using collapse(2), they become private to (i, j) iteration.
    double  tmp1, tmp2, tmp3;

    for (j = jst; j <= jend; j++) {
      // tmp1, tmp2, tmp3 are redefined multiple times within this inner loop body
      // for calculations of d, a, b, and c. Their values do not carry over
      // between (i,j) iterations, nor are they shared within this inner block
      // across the different calculation sections (d, a, b, c).

      // --- Calculations for d[i][j][0-4][0-4] ---
      // Based on u[i][j][k][...]
      tmp1 = 1.0 / u[i][j][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      d[i][j][0][0] =  1.0 + dt_2_tx1_dx1 + dt_2_ty1_dy1 + dt_2_tz1_dz1;
      d[i][j][0][1] =  0.0;
      d[i][j][0][2] =  0.0;
      d[i][j][0][3] =  0.0;
      d[i][j][0][4] =  0.0;

      d[i][j][1][0] =  dt_2 * (  tx1 * ( - r43_c34 * tmp2 * u[i][j][k][1] )
	                     + ty1 * ( -       c34 * tmp2 * u[i][j][k][1] )
	                     + tz1 * ( -       c34 * tmp2 * u[i][j][k][1] ) );
      d[i][j][1][1] =  1.0
	+ dt_2 * (  tx1 * r43_c34 * tmp1
	     + ty1 *       c34 * tmp1
	     + tz1 *       c34 * tmp1 )
	+ dt_2_tx1_dx2 + dt_2_ty1_dy2 + dt_2_tz1_dz2;
      d[i][j][1][2] = 0.0;
      d[i][j][1][3] = 0.0;
      d[i][j][1][4] = 0.0;

      d[i][j][2][0] = dt_2 * (  tx1 * ( -       c34 * tmp2 * u[i][j][k][2] )
	                     + ty1 * ( - r43_c34 * tmp2 * u[i][j][k][2] )
	                     + tz1 * ( -       c34 * tmp2 * u[i][j][k][2] ) );
      d[i][j][2][1] = 0.0;
      d[i][j][2][2] = 1.0
	+ dt_2 * (  tx1 *       c34 * tmp1
	     + ty1 * r43_c34 * tmp1
	     + tz1 *       c34 * tmp1 )
	+ dt_2_tx1_dx3 + dt_2_ty1_dy3 + dt_2_tz1_dz3;
      d[i][j][2][3] = 0.0;
      d[i][j][2][4] = 0.0;

      d[i][j][3][0] = dt_2 * (  tx1 * ( -       c34 * tmp2 * u[i][j][k][3] )
	                     + ty1 * ( -       c34 * tmp2 * u[i][j][k][3] )
	                     + tz1 * ( - r43_c34 * tmp2 * u[i][j][k][3] ) );
      d[i][j][3][1] = 0.0;
      d[i][j][3][2] = 0.0;
      d[i][j][3][3] = 1.0
	+ dt_2 * (  tx1 *       c34 * tmp1
	     + ty1 *       c34 * tmp1
	     + tz1 * r43_c34 * tmp1 )
	+ dt_2_tx1_dx4 + dt_2_ty1_dy4 + dt_2_tz1_dz4;
      d[i][j][3][4] = 0.0;

      d[i][j][4][0] = dt_2 * ( tx1 * ( - r43_c34_m_c1345 * tmp3 * pow2(u[i][j][k][1])
		                    - c34_m_c1345 * tmp3 * pow2(u[i][j][k][2])
		                    - c34_m_c1345 * tmp3 * pow2(u[i][j][k][3])
		                    - c1345 * tmp2 * u[i][j][k][4] )
	                     + ty1 * ( - c34_m_c1345 * tmp3 * pow2(u[i][j][k][1])
		                    - r43_c34_m_c1345 * tmp3 * pow2(u[i][j][k][2])
		                    - c34_m_c1345 * tmp3 * pow2(u[i][j][k][3])
		                    - c1345 * tmp2 * u[i][j][k][4] )
	                     + tz1 * ( - c34_m_c1345 * tmp3 * pow2(u[i][j][k][1])
		                    - c34_m_c1345 * tmp3 * pow2(u[i][j][k][2])
		                    - r43_c34_m_c1345 * tmp3 * pow2(u[i][j][k][3])
		                    - c1345 * tmp2 * u[i][j][k][4] ) );
      d[i][j][4][1] = dt_2 * ( tx1 * r43_c34_m_c1345 * tmp2 * u[i][j][k][1]
	                     + ty1 * c34_m_c1345 * tmp2 * u[i][j][k][1]
	                     + tz1 * c34_m_c1345 * tmp2 * u[i][j][k][1] );
      d[i][j][4][2] = dt_2 * ( tx1 * c34_m_c1345 * tmp2 * u[i][j][k][2]
	                     + ty1 * r43_c34_m_c1345 * tmp2 * u[i][j][k][2]
	                     + tz1 * c34_m_c1345 * tmp2 * u[i][j][k][2] );
      d[i][j][4][3] = dt_2 * ( tx1 * c34_m_c1345 * tmp2 * u[i][j][k][3]
	                     + ty1 * c34_m_c1345 * tmp2 * u[i][j][k][3]
	                     + tz1 * r43_c34_m_c1345 * tmp2 * u[i][j][k][3] );
      d[i][j][4][4] = 1.0
	+ dt_2 * ( tx1 * c1345 * tmp1
		       + ty1 * c1345 * tmp1
		       + tz1 * c1345 * tmp1 )
        + dt_2_tx1_dx5 + dt_2_ty1_dy5 + dt_2_tz1_dz5;


      // --- Calculations for a[i][j][0-4][0-4] ---
      // Based on u[i][j][k-1][...]
      tmp1 = 1.0 / u[i][j][k-1][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      a[i][j][0][0] = - dt_tz1 * dz1;
      a[i][j][0][1] =   0.0;
      a[i][j][0][2] =   0.0;
      a[i][j][0][3] = - dt_tz2;
      a[i][j][0][4] =   0.0;

      a[i][j][1][0] = - dt_tz2 * ( - ( u[i][j][k-1][1]*u[i][j][k-1][3] ) * tmp2 )
	- dt_tz1 * ( - c34 * tmp2 * u[i][j][k-1][1] );
      a[i][j][1][1] = - dt_tz2 * ( u[i][j][k-1][3] * tmp1 )
	- dt_tz1 * c34 * tmp1
	- dt_tz1 * dz2 ;
      a[i][j][1][2] = 0.0;
      a[i][j][1][3] = - dt_tz2 * ( u[i][j][k-1][1] * tmp1 );
      a[i][j][1][4] = 0.0;

      a[i][j][2][0] = - dt_tz2 * ( - ( u[i][j][k-1][2]*u[i][j][k-1][3] ) * tmp2 )
	- dt_tz1 * ( - c34 * tmp2 * u[i][j][k-1][2] );
      a[i][j][2][1] = 0.0;
      a[i][j][2][2] = - dt_tz2 * ( u[i][j][k-1][3] * tmp1 )
	- dt_tz1 * c34 * tmp1
	- dt_tz1 * dz3;
      a[i][j][2][3] = - dt_tz2 * ( u[i][j][k-1][2] * tmp1 );
      a[i][j][2][4] = 0.0;

      a[i][j][3][0] = - dt_tz2
	* ( - pow2( u[i][j][k-1][3] * tmp1 )
	    + c2_05
	    * ( ( pow2(u[i][j][k-1][1])
		  + pow2(u[i][j][k-1][2])
		  + pow2(u[i][j][k-1][3]) ) * tmp2 ) )
	- dt_tz1 * ( - r43_c34 * tmp2 * u[i][j][k-1][3] );
      a[i][j][3][1] = - dt_tz2 * ( - C2 * ( u[i][j][k-1][1] * tmp1 ) );
      a[i][j][3][2] = - dt_tz2 * ( - C2 * ( u[i][j][k-1][2] * tmp1 ) );
      a[i][j][3][3] = - dt_tz2 * ( 2.0 - C2 ) * ( u[i][j][k-1][3] * tmp1 )
	- dt_tz1 * r43_c34 * tmp1
	- dt_tz1 * dz4;
      a[i][j][3][4] = - dt_tz2 * C2;

      a[i][j][4][0] = - dt_tz2
	* ( ( C2 * (  pow2(u[i][j][k-1][1]) + pow2(u[i][j][k-1][2]) + pow2(u[i][j][k-1][3]) ) * tmp2
	      - C1 * ( u[i][j][k-1][4] * tmp1 ) )
	    * ( u[i][j][k-1][3] * tmp1 ) )
	- dt_tz1
	* ( - c34_m_c1345 * tmp3 * pow2(u[i][j][k-1][1])
	    - c34_m_c1345 * tmp3 * pow2(u[i][j][k-1][2])
	    - r43_c34_m_c1345 * tmp3 * pow2(u[i][j][k-1][3])
	    - c1345 * tmp2 * u[i][j][k-1][4] );
      a[i][j][4][1] = - dt_tz2 * ( - C2 * ( u[i][j][k-1][1]*u[i][j][k-1][3] ) * tmp2 )
	- dt_tz1 * c34_m_c1345 * tmp2 * u[i][j][k-1][1];
      a[i][j][4][2] = - dt_tz2 * ( - C2 * ( u[i][j][k-1][2]*u[i][j][k-1][3] ) * tmp2 )
	- dt_tz1 * c34_m_c1345 * tmp2 * u[i][j][k-1][2];
      a[i][j][4][3] = - dt_tz2 * ( C1 * ( u[i][j][k-1][4] * tmp1 )
            - c2_05 * ( (  pow2(u[i][j][k-1][1]) + pow2(u[i][j][k-1][2]) + 3.0*pow2(u[i][j][k-1][3]) ) * tmp2 ) )
	- dt_tz1 * r43_c34_m_c1345 * tmp2 * u[i][j][k-1][3];
      a[i][j][4][4] = - dt_tz2 * ( C1 * ( u[i][j][k-1][3] * tmp1 ) )
	- dt_tz1 * c1345 * tmp1
	- dt_tz1 * dz5;


      // --- Calculations for b[i][j][0-4][0-4] ---
      // Based on u[i][j-1][k][...]
      tmp1 = 1.0 / u[i][j-1][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      b[i][j][0][0] = - dt_ty1 * dy1;
      b[i][j][0][1] =   0.0;
      b[i][j][0][2] = - dt_ty2;
      b[i][j][0][3] =   0.0;
      b[i][j][0][4] =   0.0;

      b[i][j][1][0] = - dt_ty2 * ( - ( u[i][j-1][k][1]*u[i][j-1][k][2] ) * tmp2 )
	- dt_ty1 * ( - c34 * tmp2 * u[i][j-1][k][1] );
      b[i][j][1][1] = - dt_ty2 * ( u[i][j-1][k][2] * tmp1 )
	- dt_ty1 * c34 * tmp1
	- dt_ty1 * dy2;
      b[i][j][1][2] = - dt_ty2 * ( u[i][j-1][k][1] * tmp1 );
      b[i][j][1][3] = 0.0;
      b[i][j][1][4] = 0.0;

      b[i][j][2][0] = - dt_ty2
	* ( - pow2( u[i][j-1][k][2] * tmp1 )
	    + c2_05 * ( (  pow2(u[i][j-1][k][1])
			       + pow2(u[i][j-1][k][2])
			       + pow2(u[i][j-1][k][3]) )
			    * tmp2 ) )
	- dt_ty1 * ( - r43_c34 * tmp2 * u[i][j-1][k][2] );
      b[i][j][2][1] = - dt_ty2 * ( - C2 * ( u[i][j-1][k][1] * tmp1 ) );
      b[i][j][2][2] = - dt_ty2 * ( ( 2.0 - C2 ) * ( u[i][j-1][k][2] * tmp1 ) )
	- dt_ty1 * r43_c34 * tmp1
	- dt_ty1 * dy3;
      b[i][j][2][3] = - dt_ty2 * ( - C2 * ( u[i][j-1][k][3] * tmp1 ) );
      b[i][j][2][4] = - dt_ty2 * C2;

      b[i][j][3][0] = - dt_ty2 * ( - ( u[i][j-1][k][2]*u[i][j-1][k][3] ) * tmp2 )
	- dt_ty1 * ( - c34 * tmp2 * u[i][j-1][k][3] );
      b[i][j][3][1] = 0.0;
      b[i][j][3][2] = - dt_ty2 * ( u[i][j-1][k][3] * tmp1 );
      b[i][j][3][3] = - dt_ty2 * ( u[i][j-1][k][2] * tmp1 )
	- dt_ty1 * c34 * tmp1
	- dt_ty1 * dy4;
      b[i][j][3][4] = 0.0;

      b[i][j][4][0] = - dt_ty2 * ( ( C2 * (  pow2(u[i][j-1][k][1]) + pow2(u[i][j-1][k][2]) + pow2(u[i][j-1][k][3]) ) * tmp2
	      - C1 * ( u[i][j-1][k][4] * tmp1 ) )
	    * ( u[i][j-1][k][2] * tmp1 ) )
	- dt_ty1
	* ( - c34_m_c1345 * tmp3 * pow2(u[i][j-1][k][1])
	    - r43_c34_m_c1345 * tmp3 * pow2(u[i][j-1][k][2])
	    - c34_m_c1345 * tmp3 * pow2(u[i][j-1][k][3])
	    - c1345 * tmp2 * u[i][j-1][k][4] );
      b[i][j][4][1] = - dt_ty2 * ( - C2 * ( u[i][j-1][k][1]*u[i][j-1][k][2] ) * tmp2 )
	- dt_ty1 * c34_m_c1345 * tmp2 * u[i][j-1][k][1];
      b[i][j][4][2] = - dt_ty2 * ( C1 * ( u[i][j-1][k][4] * tmp1 )
	    - c2_05
	    * ( (  pow2(u[i][j-1][k][1])
                   + 3.0 * pow2(u[i][j-1][k][2])
		   + pow2(u[i][j-1][k][3]) ) * tmp2 ) )
	- dt_ty1 * r43_c34_m_c1345 * tmp2 * u[i][j-1][k][2];
      b[i][j][4][3] = - dt_ty2 * ( - C2 * ( u[i][j-1][k][2]*u[i][j-1][k][3] ) * tmp2 )
	- dt_ty1 * c34_m_c1345 * tmp2 * u[i][j-1][k][3];
      b[i][j][4][4] = - dt_ty2 * ( C1 * ( u[i][j-1][k][2] * tmp1 ) )
	- dt_ty1 * c1345 * tmp1
	- dt_ty1 * dy5;


      // --- Calculations for c[i][j][0-4][0-4] ---
      // Based on u[i-1][j][k][...]
      tmp1 = 1.0 / u[i-1][j][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      c[i][j][0][0] = - dt_tx1 * dx1;
      c[i][j][0][1] = - dt_tx2;
      c[i][j][0][2] =   0.0;
      c[i][j][0][3] =   0.0;
      c[i][j][0][4] =   0.0;

      c[i][j][1][0] = - dt_tx2
	* ( - pow2( u[i-1][j][k][1] * tmp1 )
	    + C2 * 0.50 * (  pow2(u[i-1][j][k][1]) + pow2(u[i-1][j][k][2]) + pow2(u[i-1][j][k][3]) ) * tmp2 )
	- dt_tx1 * ( - r43_c34 * tmp2 * u[i-1][j][k][1] );
      c[i][j][1][1] = - dt_tx2
	* ( ( 2.0 - C2 ) * ( u[i-1][j][k][1] * tmp1 ) )
	- dt_tx1 * r43_c34 * tmp1
	- dt_tx1 * dx2;
      c[i][j][1][2] = - dt_tx2
	* ( - C2 * ( u[i-1][j][k][2] * tmp1 ) );
      c[i][j][1][3] = - dt_tx2
	* ( - C2 * ( u[i-1][j][k][3] * tmp1 ) );
      c[i][j][1][4] = - dt_tx2 * C2;

      c[i][j][2][0] = - dt_tx2
	* ( - ( u[i-1][j][k][1] * u[i-1][j][k][2] ) * tmp2 )
	- dt_tx1 * ( - c34 * tmp2 * u[i-1][j][k][2] );
      c[i][j][2][1] = - dt_tx2 * ( u[i-1][j][k][2] * tmp1 );
      c[i][j][2][2] = - dt_tx2 * ( u[i-1][j][k][1] * tmp1 )
	- dt_tx1 * c34 * tmp1
	- dt_tx1 * dx3;
      c[i][j][2][3] = 0.0;
      c[i][j][2][4] = 0.0;

      c[i][j][3][0] = - dt_tx2
	* ( - ( u[i-1][j][k][1]*u[i-1][j][k][3] ) * tmp2 )
	- dt_tx1 * ( - c34 * tmp2 * u[i-1][j][k][3] );
      c[i][j][3][1] = - dt_tx2 * ( u[i-1][j][k][3] * tmp1 );
      c[i][j][3][2] = 0.0;
      c[i][j][3][3] = - dt_tx2 * ( u[i-1][j][k][1] * tmp1 )
	- dt_tx1 * c34 * tmp1
	- dt_tx1 * dx4;
      c[i][j][3][4] = 0.0;

      c[i][j][4][0] = - dt_tx2
	* ( ( C2 * (  pow2(u[i-1][j][k][1]) + pow2(u[i-1][j][k][2]) + pow2(u[i-1][j][k][3]) ) * tmp2
	      - C1 * ( u[i-1][j][k][4] * tmp1 ) )
	    * ( u[i-1][j][k][1] * tmp1 ) )
	- dt_tx1
	* ( - r43_c34_m_c1345 * tmp3 * pow2(u[i-1][j][k][1])
	    - c34_m_c1345 * tmp3 * pow2(u[i-1][j][k][2])
	    - c34_m_c1345 * tmp3 * pow2(u[i-1][j][k][3])
	    - c1345 * tmp2 * u[i-1][j][k][4] );
      c[i][j][4][1] = - dt_tx2 * ( C1 * ( u[i-1][j][k][4] * tmp1 )
	    - c2_05 * ( (  3.0*pow2(u[i-1][j][k][1]) + pow2(u[i-1][j][k][2]) + pow2(u[i-1][j][k][3]) ) * tmp2 ) )
	- dt_tx1 * r43_c34_m_c1345 * tmp2 * u[i-1][j][k][1];
      c[i][j][4][2] = - dt_tx2 * ( - C2 * ( u[i-1][j][k][2]*u[i-1][j][k][1] ) * tmp2 )
	- dt_tx1 * c34_m_c1345 * tmp2 * u[i-1][j][k][2];
      c[i][j][4][3] = - dt_tx2 * ( - C2 * ( u[i-1][j][k][3]*u[i-1][j][k][1] ) * tmp2 )
	- dt_tx1 * c34_m_c1345 * tmp2 * u[i-1][j][k][3];
      c[i][j][4][4] = - dt_tx2 * ( C1 * ( u[i-1][j][k][1] * tmp1 ) )
	- dt_tx1 * c1345 * tmp1
	- dt_tx1 * dx5;
    }
  }
}