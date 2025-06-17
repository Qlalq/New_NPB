static void binvcrhs(double lhs[5][5], double c[5][5], double r[5]) {
  double pivot, coeff;
  
  #pragma omp parallel
  {
    #pragma omp single
    {
      pivot = 1.00 / lhs[0][0];
      lhs[0][1] = lhs[0][1] * pivot;
      lhs[0][2] = lhs[0][2] * pivot;
      lhs[0][3] = lhs[0][3] * pivot;
      lhs[0][4] = lhs[0][4] * pivot;
      c[0][0] = c[0][0] * pivot;
      c[0][1] = c[0][1] * pivot;
      c[0][2] = c[0][2] * pivot;
      c[0][3] = c[0][3] * pivot;
      c[0][4] = c[0][4] * pivot;
      r[0] = r[0] * pivot;
    }

    #pragma omp for
    for (int i = 1; i < 5; i++) {
      coeff = lhs[i][0];
      for (int j = 1; j < 5; j++) {
        lhs[i][j] -= coeff * lhs[0][j];
      }
      for (int j = 0; j < 5; j++) {
        c[i][j] -= coeff * c[0][j];
      }
      r[i] -= coeff * r[0];
    }

    #pragma omp single
    {
      pivot = 1.00 / lhs[1][1];
      lhs[1][2] = lhs[1][2] * pivot;
      lhs[1][3] = lhs[1][3] * pivot;
      lhs[1][4] = lhs[1][4] * pivot;
      c[1][0] = c[1][0] * pivot;
      c[1][1] = c[1][1] * pivot;
      c[1][2] = c[1][2] * pivot;
      c[1][3] = c[1][3] * pivot;
      c[1][4] = c[1][4] * pivot;
      r[1] = r[1] * pivot;
    }

    #pragma omp for
    for (int i = 0; i < 5; i++) {
      if (i != 1) {
        coeff = lhs[i][1];
        for (int j = 2; j < 5; j++) {
          lhs[i][j] -= coeff * lhs[1][j];
        }
        for (int j = 0; j < 5; j++) {
          c[i][j] -= coeff * c[1][j];
        }
        r[i] -= coeff * r[1];
      }
    }

    #pragma omp single
    {
      pivot = 1.00 / lhs[2][2];
      lhs[2][3] = lhs[2][3] * pivot;
      lhs[2][4] = lhs[2][4] * pivot;
      c[2][0] = c[2][0] * pivot;
      c[2][1] = c[2][1] * pivot;
      c[2][2] = c[2][2] * pivot;
      c[2][3] = c[2][3] * pivot;
      c[2][4] = c[2][4] * pivot;
      r[2] = r[2] * pivot;
    }

    #pragma omp for
    for (int i = 0; i < 5; i++) {
      if (i != 2) {
        coeff = lhs[i][2];
        for (int j = 3; j < 5; j++) {
          lhs[i][j] -= coeff * lhs[2][j];
        }
        for (int j = 0; j < 5; j++) {
          c[i][j] -= coeff * c[2][j];
        }
        r[i] -= coeff * r[2];
      }
    }

    #pragma omp single
    {
      pivot = 1.00 / lhs[3][3];
      lhs[3][4] = lhs[3][4] * pivot;
      c[3][0] = c[3][0] * pivot;
      c[3][1] = c[3][1] * pivot;
      c[3][2] = c[3][2] * pivot;
      c[3][3] = c[3][3] * pivot;
      c[3][4] = c[3][4] * pivot;
      r[3] = r[3] * pivot;
    }

    #pragma omp for
    for (int i = 0; i < 5; i++) {
      if (i != 3) {
        coeff = lhs[i][3];
        for (int j = 4; j < 5; j++) {
          lhs[i][j] -= coeff * lhs[3][j];
        }
        for (int j = 0; j < 5; j++) {
          c[i][j] -= coeff * c[3][j];
        }
        r[i] -= coeff * r[3];
      }
    }

    #pragma omp single
    {
      pivot = 1.00 / lhs[4][4];
      c[4][0] = c[4][0] * pivot;
      c[4][1] = c[4][1] * pivot;
      c[4][2] = c[4][2] * pivot;
      c[4][3] = c[4][3] * pivot;
      c[4][4] = c[4][4] * pivot;
      r[4] = r[4] * pivot;
    }

    #pragma omp for
    for (int i = 0; i < 5; i++) {
      if (i != 4) {
        coeff = lhs[i][4];
        for (int j = 0; j < 5; j++) {
          c[i][j] -= coeff * c[4][j];
        }
        r[i] -= coeff * r[4];
      }
    }
  }
}