// Assume necessary declarations for grid_points (int array),
// u, fjac, njac, lhs (multi-dimensional double arrays),
// constants (c1, c2, c3c4, con43, c1345, dt, tx1, tx2, dx1-dx5 - double),
// AA, BB, CC (int constants/enums), and pow2 (function/macro)
// exist at file scope or are accessible in the scope where lhsx is defined.

// Forward declarations for the helper functions
static void compute_fjac_njac_plane(int j, int k);
static void compute_lhs_plane(int j, int k);


// Helper function to compute fjac and njac for a given (j, k) plane.
// This corresponds to the first inner loop (i=0 to grid_points[0]-1)
// in the original function.
static void compute_fjac_njac_plane(int j, int k) {
    // Phase 1: Compute fjac and njac for all i at this (j, k) plane.
    // This loop calculates values fjac[i][j][k] and njac[i][j][k] which depend only on u[i][j][k].
    // The computations are independent across 'i' for a fixed (j, k).
    // This inner loop COULD theoretically be parallelized itself (e.g., with #pragma omp for),
    // but parallelizing the outer j/k loops in the caller (`lhsx`) is typically more efficient
    // as it involves less overhead per work unit and avoids potential cache issues associated
    // with strided i-access across threads writing to fjac/njac for the same (j,k) plane.
    for (int i = 0; i < grid_points[0]; i++) {
        // Temporary variables used within the i-loop calculation.
        // Declaring them inside the i-loop ensures they are private to each iteration.
        // When the calling loop (j/k loop in lhsx) is parallelized, these automatics
        // within the called function are handled correctly as thread-private by OpenMP,
        // adhering to the privatization principle.
        double tmp1 = 1.0 / u[i][j][k][0];
        double tmp2 = tmp1 * tmp1;
        double tmp3 = tmp1 * tmp2;

        // Compute elements of fjac[i][j][k]
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

        // Compute elements of njac[i][j][k]
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
    } // end i loop (Phase 1: fjac and njac)
}

// Helper function to compute lhs for a given (j, k) plane.
// This corresponds to the second inner loop (i=1 to grid_points[0]-2)
// in the original function.
static void compute_lhs_plane(int j, int k) {
    // Phase 2: Compute lhs for i=1..size-2 at this (j, k) plane.
    // This loop calculates lhs[i][j][k] using stencil access (i-1, i, i+1)
    // to the fjac and njac values computed in compute_fjac_njac_plane for the same (j, k).
    // Due to the dependencies across 'i' in this phase (accessing i-1 and i+1),
    // the inner i loop is NOT directly parallelizable with a simple OpenMP #pragma omp for
    // WITHOUT specific techniques like domain decomposition with ghost zones or different algorithms.
    // However, the entire function call for a given (j, k) is independent of other (j, k) pairs,
    // making the outer j/k loops ideal for parallelization.

    // Temporary variables loop-invariant within the i-loop.
    // Declaring them here makes them local to this function call instance,
    // ensuring they are private to each thread when the caller loop (j/k) is parallelized.
    double tmp1_phase2 = dt * tx1;
    double tmp2_phase2 = dt * tx2;

    for (int i = 1; i < grid_points[0]-1; i++) {
        // Compute elements of lhs[i][j][k][AA]
        lhs[i][j][k][AA][0][0] = - tmp2_phase2 * fjac[i-1][j][k][0][0]
                                 - tmp1_phase2 * njac[i-1][j][k][0][0]
                                 - tmp1_phase2 * dx1;
        lhs[i][j][k][AA][0][1] = - tmp2_phase2 * fjac[i-1][j][k][0][1]
                                 - tmp1_phase2 * njac[i-1][j][k][0][1];
        lhs[i][j][k][AA][0][2] = - tmp2_phase2 * fjac[i-1][j][k][0][2]
                                 - tmp1_phase2 * njac[i-1][j][k][0][2];
        lhs[i][j][k][AA][0][3] = - tmp2_phase2 * fjac[i-1][j][k][0][3]
                                 - tmp1_phase2 * njac[i-1][j][k][0][3];
        lhs[i][j][k][AA][0][4] = - tmp2_phase2 * fjac[i-1][j][k][0][4]
                                 - tmp1_phase2 * njac[i-1][j][k][0][4];
        lhs[i][j][k][AA][1][0] = - tmp2_phase2 * fjac[i-1][j][k][1][0]
                                 - tmp1_phase2 * njac[i-1][j][k][1][0];
        lhs[i][j][k][AA][1][1] = - tmp2_phase2 * fjac[i-1][j][k][1][1]
                                 - tmp1_phase2 * njac[i-1][j][k][1][1]
                                 - tmp1_phase2 * dx2;
        lhs[i][j][k][AA][1][2] = - tmp2_phase2 * fjac[i-1][j][k][1][2]
                                 - tmp1_phase2 * njac[i-1][j][k][1][2];
        lhs[i][j][k][AA][1][3] = - tmp2_phase2 * fjac[i-1][j][k][1][3]
                                 - tmp1_phase2 * njac[i-1][j][k][1][3];
        lhs[i][j][k][AA][1][4] = - tmp2_phase2 * fjac[i-1][j][k][1][4]
                                 - tmp1_phase2 * njac[i-1][j][k][1][4];
        lhs[i][j][k][AA][2][0] = - tmp2_phase2 * fjac[i-1][j][k][2][0]
                                 - tmp1_phase2 * njac[i-1][j][k][2][0];
        lhs[i][j][k][AA][2][1] = - tmp2_phase2 * fjac[i-1][j][k][2][1]
                                 - tmp1_phase2 * njac[i-1][j][k][2][1];
        lhs[i][j][k][AA][2][2] = - tmp2_phase2 * fjac[i-1][j][k][2][2]
                                 - tmp1_phase2 * njac[i-1][j][k][2][2]
                                 - tmp1_phase2 * dx3;
        lhs[i][j][k][AA][2][3] = - tmp2_phase2 * fjac[i-1][j][k][2][3]
                                 - tmp1_phase2 * njac[i-1][j][k][2][3];
        lhs[i][j][k][AA][2][4] = - tmp2_phase2 * fjac[i-1][j][k][2][4]
                                 - tmp1_phase2 * njac[i-1][j][k][2][4];
        lhs[i][j][k][AA][3][0] = - tmp2_phase2 * fjac[i-1][j][k][3][0]
                                 - tmp1_phase2 * njac[i-1][j][k][3][0];
        lhs[i][j][k][AA][3][1] = - tmp2_phase2 * fjac[i-1][j][k][3][1]
                                 - tmp1_phase2 * njac[i-1][j][k][3][1];
        lhs[i][j][k][AA][3][2] = - tmp2_phase2 * fjac[i-1][j][k][3][2]
                                 - tmp1_phase2 * njac[i-1][j][k][3][2];
        lhs[i][j][k][AA][3][3] = - tmp2_phase2 * fjac[i-1][j][k][3][3]
                                 - tmp1_phase2 * njac[i-1][j][k][3][3]
                                 - tmp1_phase2 * dx4;
        lhs[i][j][k][AA][3][4] = - tmp2_phase2 * fjac[i-1][j][k][3][4]
                                 - tmp1_phase2 * njac[i-1][j][k][3][4];
        lhs[i][j][k][AA][4][0] = - tmp2_phase2 * fjac[i-1][j][k][4][0]
                                 - tmp1_phase2 * njac[i-1][j][k][4][0];
        lhs[i][j][k][AA][4][1] = - tmp2_phase2 * fjac[i-1][j][k][4][1]
                                 - tmp1_phase2 * njac[i-1][j][k][4][1];
        lhs[i][j][k][AA][4][2] = - tmp2_phase2 * fjac[i-1][j][k][4][2]
                                 - tmp1_phase2 * njac[i-1][j][k][4][2];
        lhs[i][j][k][AA][4][3] = - tmp2_phase2 * fjac[i-1][j][k][4][3]
                                 - tmp1_phase2 * njac[i-1][j][k][4][3];
        lhs[i][j][k][AA][4][4] = - tmp2_phase2 * fjac[i-1][j][k][4][4]
                                 - tmp1_phase2 * njac[i-1][j][k][4][4]
                                 - tmp1_phase2 * dx5;

        // Compute elements of lhs[i][j][k][BB]
        lhs[i][j][k][BB][0][0] = 1.0
                                + tmp1_phase2 * 2.0 * njac[i][j][k][0][0]
                                + tmp1_phase2 * 2.0 * dx1;
        lhs[i][j][k][BB][0][1] = tmp1_phase2 * 2.0 * njac[i][j][k][0][1];
        lhs[i][j][k][BB][0][2] = tmp1_phase2 * 2.0 * njac[i][j][k][0][2];
        lhs[i][j][k][BB][0][3] = tmp1_phase2 * 2.0 * njac[i][j][k][0][3];
        lhs[i][j][k][BB][0][4] = tmp1_phase2 * 2.0 * njac[i][j][k][0][4];
        lhs[i][j][k][BB][1][0] = tmp1_phase2 * 2.0 * njac[i][j][k][1][0];
        lhs[i][j][k][BB][1][1] = 1.0
                                + tmp1_phase2 * 2.0 * njac[i][j][k][1][1]
                                + tmp1_phase2 * 2.0 * dx2;
        lhs[i][j][k][BB][1][2] = tmp1_phase2 * 2.0 * njac[i][j][k][1][2];
        lhs[i][j][k][BB][1][3] = tmp1_phase2 * 2.0 * njac[i][j][k][1][3];
        lhs[i][j][k][BB][1][4] = tmp1_phase2 * 2.0 * njac[i][j][k][1][4];
        lhs[i][j][k][BB][2][0] = tmp1_phase2 * 2.0 * njac[i][j][k][2][0];
        lhs[i][j][k][BB][2][1] = tmp1_phase2 * 2.0 * njac[i][j][k][2][1];
        lhs[i][j][k][BB][2][2] = 1.0
                                + tmp1_phase2 * 2.0 * njac[i][j][k][2][2]
                                + tmp1_phase2 * 2.0 * dx3;
        lhs[i][j][k][BB][2][3] = tmp1_phase2 * 2.0 * njac[i][j][k][2][3];
        lhs[i][j][k][BB][2][4] = tmp1_phase2 * 2.0 * njac[i][j][k][2][4];
        lhs[i][j][k][BB][3][0] = tmp1_phase2 * 2.0 * njac[i][j][k][3][0];
        lhs[i][j][k][BB][3][1] = tmp1_phase2 * 2.0 * njac[i][j][k][3][1];
        lhs[i][j][k][BB][3][2] = tmp1_phase2 * 2.0 * njac[i][j][k][3][2];
        lhs[i][j][k][BB][3][3] = 1.0
                                + tmp1_phase2 * 2.0 * njac[i][j][k][3][3]
                                + tmp1_phase2 * 2.0 * dx4;
        lhs[i][j][k][BB][3][4] = tmp1_phase2 * 2.0 * njac[i][j][k][3][4];
        lhs[i][j][k][BB][4][0] = tmp1_phase2 * 2.0 * njac[i][j][k][4][0];
        lhs[i][j][k][BB][4][1] = tmp1_phase2 * 2.0 * njac[i][j][k][4][1];
        lhs[i][j][k][BB][4][2] = tmp1_phase2 * 2.0 * njac[i][j][k][4][2];
        lhs[i][j][k][BB][4][3] = tmp1_phase2 * 2.0 * njac[i][j][k][4][3];
        lhs[i][j][k][BB][4][4] = 1.0
                                + tmp1_phase2 * 2.0 * njac[i][j][k][4][4]
                                + tmp1_phase2 * 2.0 * dx5;

        // Compute elements of lhs[i][j][k][CC]
        lhs[i][j][k][CC][0][0] =  tmp2_phase2 * fjac[i+1][j][k][0][0]
                                 - tmp1_phase2 * njac[i+1][j][k][0][0]
                                 - tmp1_phase2 * dx1;
        lhs[i][j][k][CC][0][1] =  tmp2_phase2 * fjac[i+1][j][k][0][1]
                                 - tmp1_phase2 * njac[i+1][j][k][0][1];
        lhs[i][j][k][CC][0][2] =  tmp2_phase2 * fjac[i+1][j][k][0][2]
                                 - tmp1_phase2 * njac[i+1][j][k][0][2];
        lhs[i][j][k][CC][0][3] =  tmp2_phase2 * fjac[i+1][j][k][0][3]
                                 - tmp1_phase2 * njac[i+1][j][k][0][3];
        lhs[i][j][k][CC][0][4] =  tmp2_phase2 * fjac[i+1][j][k][0][4]
                                 - tmp1_phase2 * njac[i+1][j][k][0][4];
        lhs[i][j][k][CC][1][0] =  tmp2_phase2 * fjac[i+1][j][k][1][0]
                                 - tmp1_phase2 * njac[i+1][j][k][1][0];
        lhs[i][j][k][CC][1][1] =  tmp2_phase2 * fjac[i+1][j][k][1][1]
                                 - tmp1_phase2 * njac[i+1][j][k][1][1]
                                 - tmp1_phase2 * dx2;
        lhs[i][j][k][CC][1][2] =  tmp2_phase2 * fjac[i+1][j][k][1][2]
                                 - tmp1_phase2 * njac[i+1][j][k][1][2];
        lhs[i][j][k][CC][1][3] =  tmp2_phase2 * fjac[i+1][j][k][1][3]
                                 - tmp1_phase2 * njac[i+1][j][k][1][3];
        lhs[i][j][k][CC][1][4] =  tmp2_phase2 * fjac[i+1][j][k][1][4]
                                 - tmp1_phase2 * njac[i+1][j][k][1][4];
        lhs[i][j][k][CC][2][0] =  tmp2_phase2 * fjac[i+1][j][k][2][0]
                                 - tmp1_phase2 * njac[i+1][j][k][2][0];
        lhs[i][j][k][CC][2][1] =  tmp2_phase2 * fjac[i+1][j][k][2][1]
                                 - tmp1_phase2 * njac[i+1][j][k][2][1];
        lhs[i][j][k][CC][2][2] =  tmp2_phase2 * fjac[i+1][j][k][2][2]
                                 - tmp1_phase2 * njac[i+1][j][k][2][2]
                                 - tmp1_phase2 * dx3;
        lhs[i][j][k][CC][2][3] =  tmp2_phase2 * fjac[i+1][j][k][2][3]
                                 - tmp1_phase2 * njac[i+1][j][k][2][3];
        lhs[i][j][k][CC][2][4] =  tmp2_phase2 * fjac[i+1][j][k][2][4]
                                 - tmp1_phase2 * njac[i+1][j][k][2][4];
        lhs[i][j][k][CC][3][0] =  tmp2_phase2 * fjac[i+1][j][k][3][0]
                                 - tmp1_phase2 * njac[i+1][j][k][3][0];
        lhs[i][j][k][CC][3][1] =  tmp2_phase2 * fjac[i+1][j][k][3][1]
                                 - tmp1_phase2 * njac[i+1][j][k][3][1];
        lhs[i][j][k][CC][3][2] =  tmp2_phase2 * fjac[i+1][j][k][3][2]
                                 - tmp1_phase2 * njac[i+1][j][k][3][2];
        lhs[i][j][k][CC][3][3] =  tmp2_phase2 * fjac[i+1][j][k][3][3]
                                 - tmp1_phase2 * njac[i+1][j][k][3][3]
                                 - tmp1_phase2 * dx4;
        lhs[i][j][k][CC][3][4] =  tmp2_phase2 * fjac[i+1][j][k][3][4]
                                 - tmp1_phase2 * njac[i+1][j][k][3][4];
        lhs[i][j][k][CC][4][0] =  tmp2_phase2 * fjac[i+1][j][k][4][0]
                                 - tmp1_phase2 * njac[i+1][j][k][4][0];
        lhs[i][j][k][CC][4][1] =  tmp2_phase2 * fjac[i+1][j][k][4][1]
                                 - tmp1_phase2 * njac[i+1][j][k][4][1];
        lhs[i][j][k][CC][4][2] =  tmp2_phase2 * fjac[i+1][j][k][4][2]
                                 - tmp1_phase2 * njac[i+1][j][k][4][2];
        lhs[i][j][k][CC][4][3] =  tmp2_phase2 * fjac[i+1][j][k][4][3]
                                 - tmp1_phase2 * njac[i+1][j][k][4][3];
        lhs[i][j][k][CC][4][4] =  tmp2_phase2 * fjac[i+1][j][k][4][4]
                                 - tmp1_phase2 * njac[i+1][j][k][4][4]
                                 - tmp1_phase2 * dx5;
    } // end i loop (Phase 2: lhs)
}


// The main function, refactored to use helper functions for clarity and
// to highlight the independent work units (j, k planes).
static void lhsx(void) {
    int j, k;

    // Refined structure for OpenMP parallelization:
    // The outer loops over j and k define independent computation blocks (planes in the i-dimension).
    // Each (j, k) pair's computation of the i-dimension slices of
    // fjac, njac, and lhs is independent of other (j, k) pairs.
    // This structure facilitates efficient parallelization of the j and/or k loops using OpenMP.
    // By breaking the computation for each (j, k) into two distinct phases (compute fjac/njac, then compute lhs)
    // and moving them into helper functions, we encapsulate the logic for a single (j, k) plane,
    // making the parallelizable outer loops cleaner and aiding in managing temporary variables correctly.

    // Data dependencies: Dependencies exist within the i-dimension in compute_lhs_plane
    // (stencil access i-1, i, i+1). However, these dependencies are handled sequentially
    // within each call to compute_lhs_plane. There are NO loop-carried dependencies
    // across the j or k loops, meaning different iterations of j and k can safely run in parallel.

    // Load Imbalance: Assuming uniform grid size (grid_points) and relatively uniform computation
    // cost per (j, k) plane, parallelizing j and/or k provides good load balance.
    // OpenMP's scheduling clauses (like `static`, `dynamic`, `guided`) can further
    // fine-tune load distribution if needed.

    // Synchronization Overhead: Minimized by parallelizing independent j/k blocks.
    // Threads work on separate (j, k) indices, so writes to fjac, njac, and lhs
    // arrays at these indices do not conflict. Temporary variables (tmp1, tmp2, tmp3,
    // tmp1_phase2, tmp2_phase2) are local to the helper functions, ensuring thread
    // privacy and avoiding the need for synchronization like locks or atomic operations
    // for these variables.

    for (j = 1; j < grid_points[1]-1; j++) {
        for (k = 1; k < grid_points[2]-1; k++) {

            // Compute fjac and njac for the plane (j, k)
            // This function operates independently on data within the (j, k) plane.
            compute_fjac_njac_plane(j, k);

            // Compute lhs for the plane (j, k) using results from fjac and njac
            // This function operates on data within the (j, k) plane and its i-neighbors,
            // but is independent of other (j', k') planes.
            compute_lhs_plane(j, k);

        } // end k loop
    } // end j loop
}