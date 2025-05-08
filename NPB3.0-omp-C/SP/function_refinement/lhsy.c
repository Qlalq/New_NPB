// Assume necessary includes and definitions (like max, grid_points, constants, arrays) are available
// If max is not a macro, you might need to include <math.h> for fmax or define your own
// #include <stdio.h>

static void lhsy(void) {

  // Declare loop indices. Declaring within loops is good practice for scope clarity.
  // The original code declared these at the function scope.
  // int i, j, k;

  //---------------------------------------------------------------------------
  // Block 1: Calculate lhs layers 0-4 based on rho_i and vs.
  // This block updates lhs[0..4][i][j][k] based on neighbor indices (j-1, j, j+1)
  // from rho_i and vs, incorporating the intermediate calculations from the
  // original cv and rhoq arrays on the fly.
  //
  // Data Dependencies:
  // - Reads from rho_i and vs use neighbor indices (j-1, j, j+1). These are
  //   read-only dependencies on source data within an iteration, not loop-carried writes
  //   on the target array (lhs) or temporary arrays (which are now eliminated).
  // - Writes to lhs[0..4][i][j][k] are independent for different (i, j, k).
  // - The intermediate calculations equivalent to the original cv and rhoq are
  //   performed directly within the loop iteration, eliminating the loop-carried
  //   dependency on a preceding j-loop that populated temporary arrays.
  //
  // Parallelization Strategy:
  // - Fully parallelizable over i, k, and j. Can use 'collapse' clause in OpenMP
  //   across the three loops (i, k, j).
  // - Temporary scalar variables (like ru1_j_minus_1, cv_j_minus_1, etc.) are
  //   implicitly private to each thread/iteration.
  //
  // Load Imbalance:
  // - Load is balanced per (i, j, k) iteration as the work inside the inner loop
  //   is constant. The loop ranges are uniform rectangular blocks.
  //
  // Synchronization Overhead:
  // - Writes to lhs are distinct for different (i, j, k), so no synchronization
  //   (like critical sections or atomics) is needed for lhs access within this block
  //   once the temporary arrays are eliminated.
  //---------------------------------------------------------------------------
  for (int i = 1; i <= grid_points[0]-2; i++) {
    for (int k = 1; k <= grid_points[2]-2; k++) {
      // Iterate over j indices that require computation based on neighbors.
      // The range [1, grid_points[1]-2] is determined by the access patterns
      // of j-1 and j+1 relative to the original j-loop range [1, grid_points[1]-2].
      for (int j = 1; j <= grid_points[1]-2; j++) {

            // Calculate values equivalent to original cv[j-1] and rhoq[j-1]
            double ru1_j_minus_1 = c3c4 * rho_i[i][j-1][k];
            double cv_j_minus_1 = vs[i][j-1][k];
            double rhoq_j_minus_1 = max(dy3 + con43 * ru1_j_minus_1,
		      max(dy5 + c1c5*ru1_j_minus_1,
			  max(dymax + ru1_j_minus_1,
			      dy1)));

            // Calculate values equivalent to original rhoq[j]
            double ru1_j = c3c4 * rho_i[i][j][k];
            // cv_j is not directly used in the original lhs[0..4][i][j][k] update
            double rhoq_j = max(dy3 + con43 * ru1_j,
		      max(dy5 + c1c5*ru1_j,
			  max(dymax + ru1_j,
			      dy1)));

            // Calculate values equivalent to original cv[j+1] and rhoq[j+1]
            double ru1_j_plus_1 = c3c4 * rho_i[i][j+1][k];
            double cv_j_plus_1 = vs[i][j+1][k];
            double rhoq_j_plus_1 = max(dy3 + con43 * ru1_j_plus_1,
		      max(dy5 + c1c5*ru1_j_plus_1,
			  max(dymax + ru1_j_plus_1,
			      dy1)));

            // Update lhs[0..4][i][j][k] using the calculated values
	    lhs[0][i][j][k] =  0.0;
	    lhs[1][i][j][k] = -dtty2 * cv_j_minus_1 - dtty1 * rhoq_j_minus_1;
	    lhs[2][i][j][k] =  1.0 + c2dtty1 * rhoq_j;
	    lhs[3][i][j][k] =  dtty2 * cv_j_plus_1 - dtty1 * rhoq_j_plus_1;
	    lhs[4][i][j][k] =  0.0;
      }
    } // End k loop
  } // End i loop

  //---------------------------------------------------------------------------
  // Block 2: Update lhs for boundary j indices (j=1 and j=2).
  // This block specifically updates lhs layers 1-4 at j=1 and lhs layers 1-4 at j=2.
  // These updates happen *after* Block 1. Some locations updated in Block 1
  // (specifically j=1 and j=2) are further modified here.
  //
  // Data Dependencies:
  // - Writes to lhs are independent for different (i, k) pairs.
  // - Reads from lhs are the same locations being written to (read-modify-write).
  //   These are not loop-carried dependencies across (i, k).
  // - Depends on lhs values computed in Block 1 for j=1 and j=2.
  //
  // Parallelization Strategy:
  // - Parallelize the outer i and/or k loops. These iterations are independent.
  // - No specific privatization needed for scalar variables like j as it's fixed
  //   or loop indices like i, k which are handled by OpenMP.
  //
  // Load Imbalance:
  // - Load is balanced per (i, k) iteration.
  //
  // Synchronization Overhead:
  // - Writes to lhs[x][i][1][k] and lhs[x][i][2][k] are distinct for different (i, k),
  //   so no synchronization is needed for lhs access across parallel iterations
  //   within this block.
  //---------------------------------------------------------------------------
  int j_block2 = 1; // Fixed j index for operations within the i-k loops
  for (int i = 1; i <= grid_points[0]-2; i++) {
    for (int k = 1; k <= grid_points[2]-2; k++) {
      // Operations for j=1 (original j)
      lhs[2][i][j_block2][k] = lhs[2][i][j_block2][k] + comz5;
      lhs[3][i][j_block2][k] = lhs[3][i][j_block2][k] - comz4;
      lhs[4][i][j_block2][k] = lhs[4][i][j_block2][k] + comz1;
      // Operations for j=2 (original j+1)
      lhs[1][i][j_block2+1][k] = lhs[1][i][j_block2+1][k] - comz4;
      lhs[2][i][j_block2+1][k] = lhs[2][i][j_block2+1][k] + comz6;
      lhs[3][i][j_block2+1][k] = lhs[3][i][j_block2+1][k] - comz4;
      lhs[4][i][j_block2+1][k] = lhs[4][i][j_block2+1][k] + comz1;
    }
  }

  //---------------------------------------------------------------------------
  // Block 3: Update lhs for inner j indices (j from 3 to grid_points[1]-4).
  // This is a standard triply nested loop updating lhs layers 0-4.
  // The j range [3, grid_points[1]-4] does not overlap with the specific j
  // indices handled in Blocks 2 and 4, and are within the range updated by Block 1.
  //
  // Data Dependencies:
  // - Writes to lhs[0..4][i][j][k] are independent for different (i, j, k).
  // - Reads from lhs are the same locations being written to (read-modify-write).
  //   These are not loop-carried dependencies across (i, j, k).
  // - Depends on lhs values computed in Block 1 for j in [3, grid_points[1]-4].
  //
  // Parallelization Strategy:
  // - Fully parallelizable over i, j, and k. Can use 'collapse' clause in OpenMP.
  // - No specific privatization needed beyond loop indices.
  //
  // Load Imbalance:
  // - Load is balanced per (i, j, k) iteration.
  //
  // Synchronization Overhead:
  // - Writes to lhs are distinct for different (i, j, k), no synchronization needed.
  //---------------------------------------------------------------------------
  for (int i = 1; i <= grid_points[0]-2; i++) {
    for (int j_block3 = 3; j_block3 <= grid_points[1]-4; j_block3++) {
      for (int k = 1; k <= grid_points[2]-2; k++) {
	lhs[0][i][j_block3][k] = lhs[0][i][j_block3][k] + comz1;
	lhs[1][i][j_block3][k] = lhs[1][i][j_block3][k] - comz4;
	lhs[2][i][j_block3][k] = lhs[2][i][j_block3][k] + comz6;
	lhs[3][i][j_block3][k] = lhs[3][i][j_block3][k] - comz4;
	lhs[4][i][j_block3][k] = lhs[4][i][j_block3][k] + comz1;
      }
    }
  }

  //---------------------------------------------------------------------------
  // Block 4: Update lhs for boundary j indices (j=grid_points[1]-3 and j=grid_points[1]-2).
  // Similar to Block 2, handles a fixed number of j indices near the boundary.
  // The j indices here do not overlap with those in Blocks 2 and 3, but are
  // within the range updated by Block 1.
  //
  // Data Dependencies:
  // - Writes to lhs are independent for different (i, k) pairs.
  // - Reads from lhs are the same locations being written to (read-modify-write).
  //   These are not loop-carried dependencies across (i, k).
  // - Depends on lhs values computed in Block 1 for j=grid_points[1]-3 and j=grid_points[1]-2.
  //
  // Parallelization Strategy:
  // - Parallelize the outer i and/or k loops. These iterations are independent.
  // - No specific privatization needed.
  //
  // Load Imbalance:
  // - Load is balanced per (i, k) iteration.
  //
  // Synchronization Overhead:
  // - Writes to lhs[x][i][grid_points[1]-3][k] and lhs[x][i][grid_points[1]-2][k]
  //   are distinct for different (i, k), so no synchronization is needed
  //   within this block.
  //---------------------------------------------------------------------------
  int j_block4 = grid_points[1]-3; // Fixed j index for operations within the i-k loops
  for (int i = 1; i <= grid_points[0]-2; i++) {
    for (int k = 1; k <= grid_points[2]-2; k++) {
      // Operations for j=grid_points[1]-3
      lhs[0][i][j_block4][k] = lhs[0][i][j_block4][k] + comz1;
      lhs[1][i][j_block4][k] = lhs[1][i][j_block4][k] - comz4;
      lhs[2][i][j_block4][k] = lhs[2][i][j_block4][k] + comz6;
      lhs[3][i][j_block4][k] = lhs[3][i][j_block4][k] - comz4;
      // Operations for j=grid_points[1]-2 (j+1 where j=grid_points[1]-3)
      lhs[0][i][j_block4+1][k] = lhs[0][i][j_block4+1][k] + comz1;
      lhs[1][i][j_block4+1][k] = lhs[1][i][j_block4+1][k] - comz4;
      lhs[2][i][j_block4+1][k] = lhs[2][i][j_block4+1][k] + comz5;
    }
  }

  //---------------------------------------------------------------------------
  // Block 5: Calculate lhs layers 5-14 based on lhs layers 0-4 and speed.
  // This block reads from lhs layers 0-4 (which were updated in previous blocks)
  // and speed, and writes to lhs layers 5-14.
  //
  // Data Dependencies:
  // - Reads from lhs[0..4] depend on all previous blocks completing their updates
  //   to these layers. This establishes a required sequential order *between*
  //   blocks, but not within this block's loops.
  // - Reads from speed use neighbor j indices (j-1, j+1). This is a read-only
  //   dependency on source data (speed), not a loop-carried write dependency on lhs[5..14].
  // - Writes to lhs[5..14][i][j][k] are independent for different (i, j, k). These
  //   layers are only written in this block and not read by subsequent blocks shown.
  //
  // Parallelization Strategy:
  // - Fully parallelizable over i, j, and k. Can use 'collapse' clause.
  // - No specific privatization needed beyond loop indices and temporary scalars.
  //
  // Load Imbalance:
  // - Load is balanced per (i, j, k) iteration.
  //
  // Synchronization Overhead:
  // - Writes to lhs[5..14] are distinct for different (i, j, k), no synchronization
  //   needed within this block.
  //---------------------------------------------------------------------------
  for (int i = 1; i <= grid_points[0]-2; i++) {
    for (int j_block5 = 1; j_block5 <= grid_points[1]-2; j_block5++) {
      for (int k = 1; k <= grid_points[2]-2; k++) {
	lhs[0+5][i][j_block5][k]  = lhs[0][i][j_block5][k];
	lhs[1+5][i][j_block5][k]  = lhs[1][i][j_block5][k] -
	  dtty2 * speed[i][j_block5-1][k];
	lhs[2+5][i][j_block5][k]  = lhs[2][i][j_block5][k];
	lhs[3+5][i][j_block5][k]  = lhs[3][i][j_block5][k] +
	  dtty2 * speed[i][j_block5+1][k];
	lhs[4+5][i][j_block5][k] = lhs[4][i][j_block5][k];
	lhs[0+10][i][j_block5][k] = lhs[0][i][j_block5][k];
	lhs[1+10][i][j_block5][k] = lhs[1][i][j_block5][k] +
	  dtty2 * speed[i][j_block5-1][k];
	lhs[2+10][i][j_block5][k] = lhs[2][i][j_block5][k];
	lhs[3+10][i][j_block5][k] = lhs[3][i][j_block5][k] -
	  dtty2 * speed[i][j_block5+1][k];
	lhs[4+10][i][j_block5][k] = lhs[4][i][j_block5][k];
      }
    }
  }
}