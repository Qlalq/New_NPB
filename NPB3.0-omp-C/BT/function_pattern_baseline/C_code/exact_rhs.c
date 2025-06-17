static void exact_rhs(void) {
  double dtemp[5], xi, eta, zeta, dtpp;
  int m, i, j, k, ip1, im1, jp1, jm1, km1, kp1;

  #pragma omp parallel for collapse(3) private(i, j, k, m) 
  for (i = 0; i < grid_points[0]; i++) {
    for (j = 0; j < grid_points[1]; j++) {
      for (k = 0; k < grid_points[2]; k++) {
        for (m = 0; m < 5; m++) {
          forcing[i][j][k][m] = 0.0;
        }
      }
    }
  }

  #pragma omp parallel for private(i, j, k, m, xi, eta, zeta, dtpp) 
  for (j = 1; j < grid_points[1]-1; j++) {
    eta = (double)j * dnym1;
    for (k = 1; k < grid_points[2]-1; k++) {
      zeta = (double)k * dnzm1;
      for (i = 0; i < grid_points[0]; i++) {
        xi = (double)i * dnxm1;
        exact_solution(xi, eta, zeta, dtemp);
        for (m = 0; m < 5; m++) {
          ue[i][m] = dtemp[m];
        }
        dtpp = 1.0 / dtemp[0];
        for (m = 1; m <= 4; m++) {
          buf[i][m] = dtpp * dtemp[m];
        }
        cuf[i]   = buf[i][1] * buf[i][1];
        buf[i][0] = cuf[i] + buf[i][2] * buf[i][2] + buf[i][3] * buf[i][3];
        q[i] = 0.5*(buf[i][1]*ue[i][1] + buf[i][2]*ue[i][2] + buf[i][3]*ue[i][3]);
      }

      #pragma omp parallel for private(i, ip1, im1, m) 
      for (i = 1; i < grid_points[0]-1; i++) {
        im1 = i-1;
        ip1 = i+1;
        forcing[i][j][k][0] = forcing[i][j][k][0] -
          tx2*(ue[ip1][1]-ue[im1][1]) +
          dx1tx1*(ue[ip1][0]-2.0*ue[i][0]+ue[im1][0]);
        forcing[i][j][k][1] = forcing[i][j][k][1] -
          tx2 * ((ue[ip1][1]*buf[ip1][1] + c2*(ue[ip1][4]-q[ip1])) -
          (ue[im1][1]*buf[im1][1] + c2*(ue[im1][4]-q[im1]))) +
          xxcon1*(buf[ip1][1]-2.0*buf[i][1]+buf[im1][1]) +
          dx2tx1*(ue[ip1][1]-2.0*ue[i][1]+ue[im1][1]);
        forcing[i][j][k][2] = forcing[i][j][k][2] -
          tx2 * (ue[ip1][2]*buf[ip1][1]-ue[im1][2]*buf[im1][1]) +
          xxcon2*(buf[ip1][2]-2.0*buf[i][2]+buf[im1][2]) +
          dx3tx1*(ue[ip1][2]-2.0*ue[i][2]+ue[im1][2]);
        forcing[i][j][k][3] = forcing[i][j][k][3] -
          tx2*(ue[ip1][3]*buf[ip1][1]-ue[im1][3]*buf[im1][1]) +
          xxcon2*(buf[ip1][3]-2.0*buf[i][3]+buf[im1][3]) +
          dx4tx1*(ue[ip1][3]-2.0*ue[i][3]+ue[im1][3]);
        forcing[i][j][k][4] = forcing[i][j][k][4] -
          tx2*(buf[ip1][1]*(c1*ue[ip1][4]-c2*q[ip1]) -
          buf[im1][1]*(c1*ue[im1][4]-c2*q[im1])) +
          0.5*xxcon3*(buf[ip1][0]-2.0*buf[i][0]+buf[im1][0]) +
          xxcon4*(cuf[ip1]-2.0*cuf[i]+cuf[im1]) +
          xxcon5*(buf[ip1][4]-2.0*buf[i][4]+buf[im1][4]) +
          dx5tx1*(ue[ip1][4]-2.0*ue[i][4]+ue[im1][4]);
      }

      for (m = 0; m < 5; m++) {
        i = 1;
        forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
          (5.0*ue[i][m] - 4.0*ue[i+1][m] +ue[i+2][m]);
        i = 2;
        forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
          (-4.0*ue[i-1][m] + 6.0*ue[i][m] -
           4.0*ue[i+1][m] +     ue[i+2][m]);
      }
      for (m = 0; m < 5; m++) {
        #pragma omp parallel for private(i) 
        for (i = 1*3; i <= grid_points[0]-3*1-1; i++) {
          forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
            (ue[i-2][m] - 4.0*ue[i-1][m] +
             6.0*ue[i][m] - 4.0*ue[i+1][m] + ue[i+2][m]);
        }
      }

      for (m = 0; m < 5; m++) {
        i = grid_points[0]-3;
        forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
          (ue[i-2][m] - 4.0*ue[i-1][m] +
           6.0*ue[i][m] - 4.0*ue[i+1][m]);
        i = grid_points[0]-2;
        forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
          (ue[i-2][m] - 4.0*ue[i-1][m] + 5.0*ue[i][m]);
      }
    }
  }
  // Repeat similar parallelization for the other loops (for j and k)
  // Ensure proper synchronization and dependency management where needed
}