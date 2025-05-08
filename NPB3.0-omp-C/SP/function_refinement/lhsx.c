#include <stdio.h>
#include <math.h> // For fmax

// Assume these are defined globally or passed in some context
// Replace with actual definitions from the original source if available
// int grid_points[3];
// double c3c4, dx2, con43, dx5, c1c5, dxmax, dx1, dttx2, dttx1, c2dttx1;
// double comz5, comz4, comz1, comz6;
// double ****lhs; // Assuming lhs is 5D: lhs[5][gp0][gp1][gp2] or similar
// double ***rho_i, ***us, ***speed; // Assuming 3D: [gp0][gp1][gp2]

// Use a macro for max for simplicity, assumes double arguments
#define max(a, b) fmax(a, b)

// Placeholder definitions for compilation and illustration
// In a real scenario, these would come from the application's includes/globals
static int grid_points[3] = {10, 10, 10}; // Example sizes
static double c3c4=1.0, dx2=1.0, con43=1.0, dx5=1.0, c1c5=1.0, dxmax=1.0, dx1=1.0, dttx2=1.0, dttx1=1.0, c2dttx1=1.0;
static double comz5=1.0, comz4=1.0, comz1=1.0, comz6=1.0;

// Example 5D array (assuming the dimensions from context)
// In a real application, memory management would be more complex
static double lhs[15][10][10][10]; // Example size: 15 components, gp0, gp1, gp2
static double rho_i[10][10][10], us[10][10][10], speed[10][10][10]; // Example 3D arrays

// Initialize with dummy data for compilation/testing
void init_dummy_data() {
    for(int i=0; i<10; ++i)
        for(int j=0; j<10; ++j)
            for(int k=0; k<10; ++k) {
                rho_i[i][j][k] = (double)(i+j+k+1);
                us[i][j][k] = (double)(i+j+k+2);
                speed[i][j][k] = (double)(i+j+k+3);
                for(int l=0; l<15; ++l)
                    lhs[l][i][j][k] = 0.0;
            }
}


static void lhsx_refined(void) {
  int i, j, k;

  // Block 1 & 2: Compute cv, rhon and initial lhs[0..4] values
  // The computation of cv and rhon for a given (j,k) pair is independent
  // of other (j,k) pairs. The use of cv and rhon to compute lhs for a given
  // (j,k) is also independent of other (j,k) pairs.
  // Privatize cv and rhon arrays within the parallel region (e.g., per (j,k) task)
  // Parallelize over j and k.
  for (j = 1; j <= grid_points[1]-2; j++) {
    for (k = 1; k <= grid_points[2]-2; k++) {
      // Declare cv and rhon arrays for this (j,k) task.
      // Size is grid_points[0] because the i loop goes up to grid_points[0]-1
      double cv[grid_points[0]];
      double rhon[grid_points[0]];

      // Original Block 1 loop (computes cv, rhon for this (j,k))
      for (i = 0; i <= grid_points[0]-1; i++) {
	double ru1 = c3c4*rho_i[i][j][k];
	cv[i] = us[i][j][k];
	rhon[i] = max(dx2+con43*ru1,
		      max(dx5+c1c5*ru1,
			  max(dxmax+ru1,
			      dx1)));
      }

      // Original Block 2 loop (computes lhs[0..4] for this (j,k) using the computed cv, rhon)
      // Note: This loop accesses cv/rhon at i-1, i, i+1, so cv/rhon must be computed
      // for indices up to grid_points[0]-1. This is handled by the Block 1 loop range.
      for (i = 1; i <= grid_points[0]-2; i++) {
	lhs[0][i][j][k] =   0.0; // These initializations are independent
	lhs[1][i][j][k] = - dttx2 * cv[i-1] - dttx1 * rhon[i-1]; // Depends on cv/rhon at i-1
	lhs[2][i][j][k] =   1.0 + c2dttx1 * rhon[i];           // Depends on rhon at i
	lhs[3][i][j][k] =   dttx2 * cv[i+1] - dttx1 * rhon[i+1]; // Depends on cv/rhon at i+1
	lhs[4][i][j][k] =   0.0; // These initializations are independent
      }
      // cv and rhon arrays go out of scope here for this (j,k) task
    }
  }

  // Block 3: Apply specific updates for i=1 and i=2
  // These updates modify lhs values previously computed in Block 2.
  // The updates for different (j,k) pairs are independent.
  // Parallelize over j and k.
  if (1 >= 1 && 1 <= grid_points[0]-2) { // Check if i=1 is in the valid range [1, grid_points[0]-2]
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
        int i_b3_start = 1; // Original code uses i=1 for updates on i=1 and i=2

        // Updates for i = 1
        lhs[2][i_b3_start][j][k] = lhs[2][i_b3_start][j][k] + comz5;
        lhs[3][i_b3_start][j][k] = lhs[3][i_b3_start][j][k] - comz4;
        lhs[4][i_b3_start][j][k] = lhs[4][i_b3_start][j][k] + comz1;

        // Updates for i = 2 (i+1 in original code, relative to i=1)
        int i_b3_end = i_b3_start + 1; // This corresponds to i+1 in original code
         if (i_b3_end >= 1 && i_b3_end <= grid_points[0]-2) { // Check bounds for i+1 = 2
            lhs[1][i_b3_end][j][k] = lhs[1][i_b3_end][j][k] - comz4;
            lhs[2][i_b3_end][j][k] = lhs[2][i_b3_end][j][k] + comz6;
            lhs[3][i_b3_end][j][k] = lhs[3][i_b3_end][j][k] - comz4;
            lhs[4][i_b3_end][j][k] = lhs[4][i_b3_end][j][k] + comz1;
         }
      }
    }
  }


  // Block 4: Apply specific updates for i in range [3, grid_points[0]-4]
  // These updates modify lhs values. The updates for different (i,j,k) are independent.
  // Parallelize over i, j, and k (or a combination).
  // Swapping loop order to parallelize j and k first is consistent with other blocks.
  for (j = 1; j <= grid_points[1]-2; j++) {
    for (k = 1; k <= grid_points[2]-2; k++) {
      for (i = 3; i <= grid_points[0]-4; i++) {
        // Check if i is in the valid range [1, grid_points[0]-2] for lhs writes
        // The original loop bounds [3, grid_points[0]-4] are within [1, grid_points[0]-2]
        // if grid_points[0] >= 5. Assuming this condition holds.
        lhs[0][i][j][k] = lhs[0][i][j][k] + comz1;
        lhs[1][i][j][k] = lhs[1][i][j][k] - comz4;
        lhs[2][i][j][k] = lhs[2][i][j][k] + comz6;
        lhs[3][i][j][k] = lhs[3][i][j][k] - comz4;
        lhs[4][i][j][k] = lhs[4][i][j][k] + comz1;
      }
    }
  }

  // Block 5: Apply specific updates for i=grid_points[0]-3 and i=grid_points[0]-2
  // These updates modify lhs values. The updates for different (j,k) are independent.
  // Parallelize over j and k.
  // Check if grid_points[0]-3 is in the valid range [1, grid_points[0]-2] for lhs writes
  if (grid_points[0]-3 >= 1 && grid_points[0]-3 <= grid_points[0]-2) {
     for (j = 1; j <= grid_points[1]-2; j++) {
        for (k = 1; k <= grid_points[2]-2; k++) {
           int i_b5_start = grid_points[0]-3; // Original code uses i = grid_points[0]-3

           // Updates for i = grid_points[0]-3
           lhs[0][i_b5_start][j][k] = lhs[0][i_b5_start][j][k] + comz1;
           lhs[1][i_b5_start][j][k] = lhs[1][i_b5_start][j][k] - comz4;
           lhs[2][i_b5_start][j][k] = lhs[2][i_b5_start][j][k] + comz6;
           lhs[3][i_b5_start][j][k] = lhs[3][i_b5_start][j][k] - comz4;

           // Updates for i = grid_points[0]-2 (i+1 in original code relative to grid_points[0]-3)
           int i_b5_end = i_b5_start + 1; // This corresponds to i+1 in original code
            if (i_b5_end >= 1 && i_b5_end <= grid_points[0]-2) { // Check bounds for i+1 = grid_points[0]-2
               lhs[0][i_b5_end][j][k] = lhs[0][i_b5_end][j][k] + comz1;
               lhs[1][i_b5_end][j][k] = lhs[1][i_b5_end][j][k] - comz4;
               lhs[2][i_b5_end][j][k] = lhs[2][i_b5_end][j][k] + comz5;
            }
        }
     }
  }


  // Block 6: Compute lhs[5..14] values
  // These computations use the final lhs[0..4] values (modified by Blocks 3,4,5) and speed.
  // The computations for different (i,j,k) are independent.
  // Parallelize over i, j, and k (or a combination).
  // Swapping loop order to parallelize j and k first is consistent with other blocks.
  for (j = 1; j <= grid_points[1]-2; j++) {
    for (k = 1; k <= grid_points[2]-2; k++) {
      for (i = 1; i <= grid_points[0]-2; i++) {
        // Check if i-1 and i+1 indices for speed and lhs are within bounds.
        // The loop for i is [1, grid_points[0]-2], so i-1 is [0, grid_points[0]-3]
        // and i+1 is [2, grid_points[0]-1]. These access ranges are assumed valid for
        // the speed and lhs arrays based on context, typically speed/lhs would be
        // sized at least grid_points[0].
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

/*
// Example of how to call the refined function (if not static)
int main() {
    init_dummy_data(); // Initialize placeholder data
    lhsx_refined();
    // You can add verification code here if needed
    printf("Refined lhsx function executed.\n");
    return 0;
}
*/