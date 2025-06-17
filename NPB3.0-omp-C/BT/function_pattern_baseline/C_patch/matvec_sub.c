static void matvec_sub(double ablock[5][5], double avec[5], double bvec[5]) {
  int i;
#pragma omp parallel for
  for (i = 0; i < 5; i++) {
    bvec[i] = bvec[i] - ablock[i][0]*avec[0]
      - ablock[i][1]*avec[1]
      - ablock[i][2]*avec[2]
      - ablock[i][3]*avec[3]
      - ablock[i][4]*avec[4];
  }
}