static void adi(void) {
    // This function represents the sequential phases of an ADI method.
    // The primary data dependencies are between these sequential phases:
    // compute_rhs -> x_solve -> y_solve -> z_solve -> add.
    // These stage-level dependencies require sequential execution of the blocks below.
    // Parallelization with OpenMP will primarily be applied *within* each of these function calls,
    // typically across independent dimensions (planes or grid points).

    // Phase 1: Compute Right-Hand Side (RHS)
    // Potential for data-parallelism over grid points.
    // Minimal synchronization needed if computations are point-wise independent.
    // Load balancing depends on the loop structure and scheduling within compute_rhs.
    {
        compute_rhs();
    }

    // Phase 2: X-Direction Solve
    // Solves independent 1D systems along X for each YZ plane.
    // Data dependency: Requires RHS computed in Phase 1. Output is input for Phase 3.
    // Potential for parallelization over YZ planes (looping over j and k indices).
    // Each 1D tridiagonal solve within this loop is typically sequential (has loop-carried dependencies).
    // Load balancing depends on the distribution of YZ planes across threads.
    // Synchronization is needed to ensure Phase 1 is complete before starting, and Phase 2 before Phase 3.
    {
        x_solve();
    }

    // Phase 3: Y-Direction Solve
    // Solves independent 1D systems along Y for each XZ plane.
    // Data dependency: Requires results from Phase 2. Output is input for Phase 4.
    // Potential for parallelization over XZ planes (looping over i and k indices).
    // Each 1D solve is typically sequential.
    // Load balancing depends on the distribution of XZ planes.
    // Synchronization: Needs Phase 2 complete.
    {
        y_solve();
    }

    // Phase 4: Z-Direction Solve
    // Solves independent 1D systems along Z for each XY plane.
    // Data dependency: Requires results from Phase 3. Output is input for Phase 5.
    // Potential for parallelization over XY planes (looping over i and j indices).
    // Each 1D solve is typically sequential.
    // Load balancing depends on the distribution of XY planes.
    // Synchronization: Needs Phase 3 complete.
    {
        z_solve();
    }

    // Phase 5: Update Solution / Add Increment
    // Data dependency: Requires results from Phase 4.
    // Potential for data-parallelism over grid points.
    // Synchronization might be needed if updating shared state beyond the main grid, or for reductions.
    // Load balancing depends on the loop structure within add.
    // Synchronization: Needs Phase 4 complete.
    {
        add();
    }
}