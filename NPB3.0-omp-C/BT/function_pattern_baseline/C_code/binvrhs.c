static void binvrhs(double lhs[5][5], double r[5]) {
  double pivot, coeff;
  
  // First pivot
  pivot = 1.00 / lhs[0][0];
  #pragma omp parallel for
  for (int i = 1; i < 5; i++) {
    lhs[0][i] = lhs[0][i] * pivot;
  }
  r[0] = r[0] * pivot;

  #pragma omp parallel for
  for (int i = 1; i < 5; i++) {
    coeff = lhs[i][0];
    for (int j = 1; j < 5; j++) {
      lhs[i][j] -= coeff * lhs[0][j];
    }
    r[i] -= coeff * r[0];
  }

  // Second pivot
  pivot = 1.00 / lhs[1][1];
  #pragma omp parallel for
  for (int i = 2; i < 5; i++) {
    lhs[1][i] = lhs[1][i] * pivot;
  }
  r[1] = r[1] * pivot;

  #pragma omp parallel for
  for (int i = 0; i < 5; i++) {
    coeff = lhs[i][1];
    for (int j = 2; j < 5; j++) {
      lhs[i][j] -= coeff * lhs[1][j];
    }
    r[i] -= coeff * r[1];
  }

  // Third pivot
  pivot = 1.00 / lhs[2][2];
  #pragma omp parallel for
  for (int i = 3; i < 5; i++) {
    lhs[2][i] = lhs[2][i] * pivot;
  }
  r[2] = r[2] * pivot;

  #pragma omp parallel for
  for (int i = 0; i < 5; i++) {
    coeff = lhs[i][2];
    for (int j = 3; j < 5; j++) {
      lhs[i][j] -= coeff * lhs[2][j];
    }
    r[i] -= coeff * r[2];
  }

  // Fourth pivot
  pivot = 1.00 / lhs[3][3];
  lhs[3][4] = lhs[3][4] * pivot;
  r[3] = r[3] * pivot;

  #pragma omp parallel for
  for (int i = 0; i < 5; i++) {
    coeff = lhs[i][3];
    lhs[i][4] -= coeff * lhs[3][4];
    r[i] -= coeff * r[3];
  }

  // Fifth pivot
  pivot = 1.00 / lhs[4][4];
  r[4] = r[4] * pivot;

  #pragma omp parallel for
  for (int i = 0; i < 4; i++) {
    coeff = lhs[i][4];
    r[i] -= coeff * r[4];
  }
}