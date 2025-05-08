static void pintgr(void) {
  int i, j, k;
  int ibeg, ifin, ifin1;
  int jbeg, jfin, jfin1;
  int iglob, iglob1, iglob2;
  int jglob, jglob1, jglob2;
  // phi1 and phi2 are reused for different purposes and index meanings
  // in sequential blocks, making them temporary for each block.
  // Declaration scope implies they persist across blocks, but re-initialization
  // effectively makes them temporary for each phase of calculation (frc1, frc2, frc3).
  // Their size [ISIZ2+2][ISIZ3+2] must accommodate the maximum ranges of (i,j), (i,k), (j,k)
  // used in the calculations.
  double phi1[ISIZ2+2][ISIZ3+2];
  double phi2[ISIZ2+2][ISIZ3+2];
  double frc1, frc2, frc3;
  // Assume pow2 is defined elsewhere, e.g., #define pow2(x) ((x)*(x))

  // 1. Calculate boundaries (sequential setup)
  // These calculations determine the loop ranges and are independent of the main loops.
  ibeg = nx;
  ifin = 0;
  iglob1 = -1;
  iglob2 = nx-1;
  if (iglob1 >= ii1 && iglob2 < ii2+nx) ibeg = 0;
  if (iglob1 >= ii1-nx && iglob2 <= ii2) ifin = nx;
  if (ii1 >= iglob1 && ii1 <= iglob2) ibeg = ii1;
  if (ii2 >= iglob1 && ii2 <= iglob2) ifin = ii2;

  jbeg = ny;
  jfin = -1;
  jglob1 = 0;
  jglob2 = ny-1;
  if (jglob1 >= ji1 && jglob2 < ji2+ny) jbeg = 0;
  if (jglob1 > ji1-ny && jglob2 <= ji2) jfin = ny;
  if (ji1 >= jglob1 && ji1 <= jglob2) jbeg = ji1;
  if (ji2 >= jglob1 && ji2 <= jglob2) jfin = ji2;

  ifin1 = ifin;
  jfin1 = jfin;
  if (ifin1 == ii2) ifin1 = ifin -1;
  if (jfin1 == ji2) jfin1 = jfin -1;

  // iglob, jglob are loop-private variables used for checks inside loops,
  // not affecting dependencies outside their iteration.

  // 2. Calculate frc1 (Integral over i-j plane)
  // This block is largely independent from the frc2 and frc3 blocks,
  // allowing for parallelization within this block's loops.

  // 2a. Initialize phi1 and phi2 for frc1 calculation
  // This loop is embarrassingly parallel. Each iteration writes to distinct elements.
  for (i = 0; i <= ISIZ2+1; i++) {
    for (k = 0; k <= ISIZ3+1; k++) {
      phi1[i][k] = 0.0; // Note: k is used here, different from the computation loop
      phi2[i][k] = 0.0; // where j is used for the second index. Assumes bounds alignment.
    }
  }

  // 2b. Compute phi1[i][j] and phi2[i][j] based on u[i][j][k] for frc1
  // This loop is parallel over i and j. Each (i, j) iteration writes to
  // a distinct phi1[i][j] and phi2[i][j] element. Reads from u are independent.
  for (i = ibeg; i <= ifin; i++) {
    iglob = i; // Loop-private variable
    for (j = jbeg; j <= jfin; j++) {
      jglob = j; // Loop-private variable
      k = ki1; // ki1 and ki2 define the range in k
      phi1[i][j] = C2*(  u[i][j][k][4]
            - 0.50 * (  pow2(u[i][j][k][1])
                  + pow2(u[i][j][k][2])
                  + pow2(u[i][j][k][3]) )
            / u[i][j][k][0] );
      k = ki2;
      phi2[i][j] = C2*(  u[i][j][k][4]
            - 0.50 * (  pow2(u[i][j][k][1])
                  + pow2(u[i][j][k][2])
                  + pow2(u[i][j][k][3]) )
            / u[i][j][k][0] );
    }
  }
  // Implicit barrier: Computation of phi1/phi2 must complete before using them in frc1 reduction.

  // 2c. Reduce phi values to compute frc1
  // This loop is a reduction operation on frc1.
  // Each (i, j) iteration computes a value independently before adding to frc1.
  // Can be parallelized using OpenMP reduction(+:frc1).
  frc1 = 0.0;
  for (i = ibeg; i <= ifin1; i++) {
    for (j = jbeg; j <= jfin1; j++) {
      frc1 = frc1 + (  phi1[i][j]
             + phi1[i+1][j] // Reads from elements computed in 2b
             + phi1[i][j+1] // Reads from elements computed in 2b
             + phi1[i+1][j+1] // Reads from elements computed in 2b
             + phi2[i][j]
             + phi2[i+1][j]
             + phi2[i][j+1]
             + phi2[i+1][j+1] );
    }
  }
  // This scalar multiplication is sequential.
  frc1 = dxi * deta * frc1;


  // 3. Calculate frc2 (Integral over i-k plane at jbeg/jfin boundaries)
  // This block is independent from the frc1 and frc3 blocks.

  // 3a. Re-initialize phi1 and phi2 for frc2 calculation
  // This loop is embarrassingly parallel.
  for (i = 0; i <= ISIZ2+1; i++) {
    for (k = 0; k <= ISIZ3+1; k++) {
      phi1[i][k] = 0.0;
      phi2[i][k] = 0.0;
    }
  }

  // 3b. Compute phi1[i][k] and phi2[i][k] based on u[i][jbeg/jfin][k] for frc2
  // These loops are parallel over i and k. Each (i, k) iteration writes to
  // a distinct phi1[i][k] or phi2[i][k] element. Reads from u are independent.
  jglob = jbeg; // Use jglob for conditional check, loop uses jbeg directly
  if (jglob == ji1) { // Check if jbeg is at the ji1 boundary
    for (i = ibeg; i <= ifin; i++) {
      iglob = i; // Loop-private variable
      for (k = ki1; k <= ki2; k++) {
  phi1[i][k] = C2*(  u[i][jbeg][k][4]
            - 0.50 * (  pow2(u[i][jbeg][k][1])
                  + pow2(u[i][jbeg][k][2])
                  + pow2(u[i][jbeg][k][3]) )
            / u[i][jbeg][k][0] );
      }
    }
  }
  jglob = jfin; // Use jglob for conditional check, loop uses jfin directly
  if (jglob == ji2) { // Check if jfin is at the ji2 boundary
    for (i = ibeg; i <= ifin; i++) {
      iglob = i; // Loop-private variable
      for (k = ki1; k <= ki2; k++) {
  phi2[i][k] = C2*(  u[i][jfin][k][4] // Note: Accessing u at jfin
            - 0.50 * (  pow2(u[i][jfin][k][1])
                  + pow2(u[i][jfin][k][2])
                  + pow2(u[i][jfin][k][3]) )
            / u[i][jfin][k][0] );
      }
    }
  }
  // Implicit barrier: Computation of phi1/phi2 must complete before using them in frc2 reduction.

  // 3c. Reduce phi values to compute frc2
  // This loop is a reduction operation on frc2.
  // Each (i, k) iteration computes a value independently before adding to frc2.
  // Can be parallelized using OpenMP reduction(+:frc2).
  frc2 = 0.0;
  for (i = ibeg; i <= ifin1; i++) {
    for (k = ki1; k <= ki2-1; k++) {
      frc2 = frc2 + (  phi1[i][k] // Reads from elements computed in 3b
             + phi1[i+1][k] // Reads from elements computed in 3b
             + phi1[i][k+1] // Reads from elements computed in 3b
             + phi1[i+1][k+1] // Reads from elements computed in 3b
             + phi2[i][k]
             + phi2[i+1][k]
             + phi2[i][k+1]
             + phi2[i+1][k+1] );
    }
  }
  // This scalar multiplication is sequential.
  frc2 = dxi * dzeta * frc2;


  // 4. Calculate frc3 (Integral over j-k plane at ibeg/ifin boundaries)
  // This block is independent from the frc1 and frc2 blocks.

  // 4a. Re-initialize phi1 and phi2 for frc3 calculation
  // This loop is embarrassingly parallel.
  for (i = 0; i <= ISIZ2+1; i++) { // Note: i and k loop bounds used for re-init
    for (k = 0; k <= ISIZ3+1; k++) {
      phi1[i][k] = 0.0;
      phi2[i][k] = 0.0;
    }
  }

  // 4b. Compute phi1[j][k] and phi2[j][k] based on u[ibeg/ifin][j][k] for frc3
  // These loops are parallel over j and k. Each (j, k) iteration writes to
  // a distinct phi1[j][k] or phi2[j][k] element. Reads from u are independent.
  iglob = ibeg; // Use iglob for conditional check, loop uses ibeg directly
  if (iglob == ii1) { // Check if ibeg is at the ii1 boundary
    // Note: phi is indexed as phi1[j][k], assuming j range maps to first dim, k to second.
    // This reuses the phi array with different index mapping assumption than before.
    for (j = jbeg; j <= jfin; j++) {
      jglob = j; // Loop-private variable
      for (k = ki1; k <= ki2; k++) {
  phi1[j][k] = C2*(  u[ibeg][j][k][4] // Note: Accessing u at ibeg
            - 0.50 * (  pow2(u[ibeg][j][k][1])
                  + pow2(u[ibeg][j][k][2])
                  + pow2(u[ibeg][j][k][3]) )
            / u[ibeg][j][k][0] );
      }
    }
  }
  iglob = ifin; // Use iglob for conditional check, loop uses ifin directly
  if (iglob == ii2) { // Check if ifin is at the ii2 boundary
    // Note: phi is indexed as phi2[j][k], assuming j range maps to first dim, k to second.
    for (j = jbeg; j <= jfin; j++) {
      jglob = j; // Loop-private variable
      for (k = ki1; k <= ki2; k++) {
  phi2[j][k] = C2*(  u[ifin][j][k][4] // Note: Accessing u at ifin
            - 0.50 * (  pow2(u[ifin][j][k][1])
                  + pow2(u[ifin][j][k][2])
                  + pow2(u[ifin][j][k][3]) )
            / u[ifin][j][k][0] );
      }
    }
  }
  // Implicit barrier: Computation of phi1/phi2 must complete before using them in frc3 reduction.

  // 4c. Reduce phi values to compute frc3
  // This loop is a reduction operation on frc3.
  // Each (j, k) iteration computes a value independently before adding to frc3.
  // Can be parallelized using OpenMP reduction(+:frc3).
  // Note: phi is indexed as phi[j][k] here.
  frc3 = 0.0;
  for (j = jbeg; j <= jfin1; j++) {
    for (k = ki1; k <= ki2-1; k++) {
      frc3 = frc3 + (  phi1[j][k] // Reads from elements computed in 4b
             + phi1[j+1][k] // Reads from elements computed in 4b
             + phi1[j][k+1] // Reads from elements computed in 4b
             + phi1[j+1][k+1] // Reads from elements computed in 4b
             + phi2[j][k]
             + phi2[j+1][k]
             + phi2[j][k+1]
             + phi2[j+1][k+1] );
    }
  }
  // This scalar multiplication is sequential.
  frc3 = deta * dzeta * frc3;

  // 5. Final summation (sequential)
  // Depends on frc1, frc2, frc3 being fully computed.
  frc = 0.25 * ( frc1 + frc2 + frc3 );
}