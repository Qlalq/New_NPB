static void tzetar(void) {
  int i, j, k;
  double t1, t2, t3;
  double ac, xvel, yvel, zvel, r1, r2, r3, r4, r5;
  double acinv, ac2u, uzik1, btuz;

  // Analysis for OpenMP parallelization:
  // 1. Data Dependencies:
  //    The calculations for rhs[m][i][j][k] within the innermost loop depend only
  //    on input values and previously calculated temporary variables (t1, t2, t3)
  //    which are based on values read at the *same* indices (i, j, k).
  //    There are no loop-carried dependencies across the i, j, or k loops.
  //    Each iteration (i, j, k) is computationally independent of other iterations.
  //
  // 2. Load Imbalance:
  //    The computation performed inside the innermost loop is largely uniform
  //    across all (i, j, k) tuples within the specified bounds.
  //    A standard parallel loop construct (e.g., #pragma omp parallel for)
  //    over any of the loops (i, j, or k) or a collapsed loop over multiple
  //    dimensions (#pragma omp parallel for collapse(N)) would distribute
  //    the work reasonably evenly, especially on a regular grid.
  //
  // 3. Synchronization Overhead:
  //    Since each iteration (i, j, k) writes to a unique set of locations
  //    in the 'rhs' array that are not accessed by other iterations, no
  //    explicit synchronization (like locks, atomics, or critical sections)
  //    is needed for the writes to 'rhs'.
  //    The temporary variables (xvel, yvel, ..., t3, etc.) are computed and
  //    used entirely within a single iteration and do not need synchronization;
  //    they can be made private to each thread/iteration by OpenMP directives.
  //
  // Conclusion:
  // The current structure of the nested loops with element-wise computations
  // is already well-suited for efficient OpenMP parallelization.
  // The loops can be parallelized directly without significant structural changes
  // to break dependencies or balance load, as the iterations are independent
  // and the workload per iteration is uniform.

  for (i = 1; i <= grid_points[0]-2; i++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (k = 1; k <= grid_points[2]-2; k++) {

        // Variables declared inside the loop are naturally private to each iteration.
        // Variables declared outside but used only for calculation within the loop
        // iteration (like the ones declared at the function scope) can be made
        // private when using OpenMP directives (e.g., using the 'private' clause).

        // Read data for current cell (i, j, k)
        xvel = us[i][j][k];
        yvel = vs[i][j][k];
        zvel = ws[i][j][k];
        ac   = speed[i][j][k];
        acinv = ainv[i][j][k];
        ac2u = ac * ac;

        r1 = rhs[0][i][j][k];
        r2 = rhs[1][i][j][k];
        r3 = rhs[2][i][j][k];
        r4 = rhs[3][i][j][k];
        r5 = rhs[4][i][j][k];

        uzik1 = u[0][i][j][k];

        // Perform calculations for current cell
        btuz  = bt * uzik1;

        t1 = btuz * acinv * (r4 + r5);
        t2 = r3 + t1;
        t3 = btuz * (r4 - r5);

        // Write results back to rhs for current cell (i, j, k)
        // These writes are to distinct memory locations for each unique (i,j,k) tuple
        // processed by the loops.
        rhs[0][i][j][k] = t2;
        rhs[1][i][j][k] = -uzik1 * r2 + xvel * t2;
        rhs[2][i][j][k] =  uzik1 * r1 + yvel * t2;
        rhs[3][i][j][k] =  zvel * t2  + t3;
        rhs[4][i][j][k] =  uzik1 * (-xvel * r2 + yvel * r1) +
          qs[i][j][k] * t2 + c2iv * ac2u * t1 + zvel * t3;
      }
    }
  }
}