static void exact_solution(double xi, double eta, double zeta,
			   double dtemp[5]) {
  int m;
  // Analysis for OpenMP refinement:
  // 1. Data Dependencies: The calculation for each dtemp[m] depends only
  //    on the inputs xi, eta, zeta, and values from the m-th column of the
  //    global/external array 'ce'. There are no loop-carried dependencies
  //    where the result of one iteration (m) depends on the result of a previous
  //    iteration (m-1) or a future iteration (m+1). Each iteration writes to
  //    a unique element of the 'dtemp' array.
  //    The loop iterations are inherently independent. No decoupling transformation is needed.
  //
  // 2. Load Imbalance: The computation performed inside the loop is the same
  //    for each value of 'm'. The number of floating-point operations is constant
  //    per iteration. This loop has a perfectly balanced workload across iterations.
  //    No structural changes for load balancing are needed.
  //
  // 3. Synchronization Overhead: Each iteration writes to a distinct memory location
  //    in the 'dtemp' array (dtemp[0], dtemp[1], ..., dtemp[4]). Assuming 'ce' is
  //    read-only within this function, there are no shared data writes between
  //    iterations that would require synchronization (like locks, atomics, or critical sections).
  //    The only synchronization needed would be the implicit barrier at the end of
  //    an OpenMP parallel loop, which is standard.
  //
  // Conclusion: The current loop structure is already optimally structured for
  // parallel execution using OpenMP. A simple #pragma omp for directive applied
  // to the loop is sufficient and efficient, as the iterations are independent,
  // balanced, and write to distinct memory locations.
  // No further structural modification is required based on the stated refinement principles.

  for (m = 0; m < 5; m++) {
    dtemp[m] =  ce[0][m] +
      xi*(ce[1][m] + xi*(ce[4][m] +
			     xi*(ce[7][m] + xi*ce[10][m]))) +
      eta*(ce[2][m] + eta*(ce[5][m] +
			       eta*(ce[8][m] + eta*ce[11][m])))+
      zeta*(ce[3][m] + zeta*(ce[6][m] +
				 zeta*(ce[9][m] +
				       zeta*ce[12][m])));
  }
}