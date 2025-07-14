#include "npb-C.h"
#include "header.h"
static void add(void);
static void adi(void);
static void error_norm(double rms[5]);
static void rhs_norm(double rms[5]);
static void exact_rhs(void);
static void exact_solution(double xi, double eta, double zeta,
			   double dtemp[5]);
static void initialize(void);
static void lhsinit(void);
static void lhsx(void);
static void lhsy(void);
static void lhsz(void);
static void compute_rhs(void);
static void set_constants(void);
static void verify(int no_time_steps, char *class, boolean *verified);
static void x_solve(void);
static void x_backsubstitute(void);
static void x_solve_cell(void);
static void matvec_sub(double ablock[5][5], double avec[5], double bvec[5]);
static void matmul_sub(double ablock[5][5], double bblock[5][5],
		       double cblock[5][5]);
static void binvcrhs(double lhs[5][5], double c[5][5], double r[5]);
static void binvrhs(double lhs[5][5], double r[5]);
static void y_solve(void);
static void y_backsubstitute(void);
static void y_solve_cell(void);
static void z_solve(void);
static void z_backsubstitute(void);
static void z_solve_cell(void);
int main(int argc, char **argv) {
  int niter, step, n3;
  int nthreads = 1;
  double navg, mflops;
  double tmax;
  boolean verified;
  char class;
  FILE *fp;
  printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version"
	 " - BT Benchmark\n\n");
  fp = fopen("inputbt.data", "r");
  if (fp != NULL) {
    printf(" Reading from input file inputbt.data");
    fscanf(fp, "%d", &niter);
    while (fgetc(fp) != '\n');
    fscanf(fp, "%lg", &dt);
    while (fgetc(fp) != '\n');
    fscanf(fp, "%d%d%d",
	   &grid_points[0],  &grid_points[1],  &grid_points[2]);
    fclose(fp);
  } else {
    printf(" No input file inputbt.data. Using compiled defaults\n");
    niter = NITER_DEFAULT;
    dt    = DT_DEFAULT;
    grid_points[0] = PROBLEM_SIZE;
    grid_points[1] = PROBLEM_SIZE;
    grid_points[2] = PROBLEM_SIZE;
  }
  printf(" Size: %3dx%3dx%3d\n",
	 grid_points[0], grid_points[1], grid_points[2]);
  printf(" Iterations: %3d   dt: %10.6f\n", niter, dt);
  if (grid_points[0] > IMAX ||
      grid_points[1] > JMAX ||
      grid_points[2] > KMAX) {
    printf(" %dx%dx%d\n", grid_points[0], grid_points[1], grid_points[2]);
    printf(" Problem size too big for compiled array sizes\n");
    exit(1);
  }

    set_constants();
    initialize();
    lhsinit();
    exact_rhs();
    adi();
    initialize();

  timer_clear(1);
  timer_start(1);



{
    for (step = 1; step <= niter; step++) {
            if (step % 20 == 0 || step == 1) {
                printf(" Time step %4d\n", step);
            }
        adi();
    }
}


  {   
    #if defined(_OPENMP)
    nthreads = omp_get_num_threads();
    #endif 
  } 

  timer_stop(1);
  tmax = timer_read(1);

  {
    verify(niter, &class, &verified);
  }

  n3 = grid_points[0]*grid_points[1]*grid_points[2];
  navg = (grid_points[0]+grid_points[1]+grid_points[2])/3.0;
  if ( tmax != 0.0 ) {
    mflops = 1.0e-6*(double)niter*
	(3478.8*(double)n3-17655.7*pow2(navg)+28023.7*navg) / tmax;
  } else {
    mflops = 0.0;
  }

  c_print_results("BT", class, grid_points[0], 
		  grid_points[1], grid_points[2], niter, nthreads,
		  tmax, mflops, "          floating point", 
		  verified, NPBVERSION,COMPILETIME, CS1, CS2, CS3, CS4, CS5, 
		  CS6, "(none)");
}
static void add(void) {
  int i, j, k, m;
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	for (m = 0; m < 5; m++) {
	  u[i][j][k][m] = u[i][j][k][m] + rhs[i][j][k][m];
	}
      }
    }
  }
}
static void adi(void) {
    compute_rhs();
    x_solve();
    y_solve();
    z_solve();
    add();
}
static void error_norm(double rms[5]) {
  int i, j, k, m, d;
  double xi, eta, zeta, u_exact[5], add;
  for (m = 0; m < 5; m++) {
    rms[m] = 0.0;
  }
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      for (k = 0; k < grid_points[2]; k++) {
	zeta = (double)k * dnzm1;
	exact_solution(xi, eta, zeta, u_exact);
	for (m = 0; m < 5; m++) {
	  add = u[i][j][k][m] - u_exact[m];
	  rms[m] = rms[m] + add*add;
	}
      }
    }
  }
  for (m = 0; m < 5; m++) {
    for (d = 0; d <= 2; d++) {
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    }
    rms[m] = sqrt(rms[m]);
  }
}
static void rhs_norm(double rms[5]) {
  int i, j, k, d, m;
  double add;
  for (m = 0; m < 5; m++) {
    rms[m] = 0.0;
  }
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	for (m = 0; m < 5; m++) {
	  add = rhs[i][j][k][m];
	  rms[m] = rms[m] + add*add;
	}
      }
    }
  }
  for (m = 0; m < 5; m++) {
    for (d = 0; d <= 2; d++) {
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    }
    rms[m] = sqrt(rms[m]);
  }
}
static void exact_rhs(void) {
{
  double dtemp[5], xi, eta, zeta, dtpp;
  int m, i, j, k, ip1, im1, jp1, jm1, km1, kp1;
  for (i = 0; i < grid_points[0]; i++) {
    for (j = 0; j < grid_points[1]; j++) {
      for (k = 0; k < grid_points[2]; k++) {
	for (m = 0; m < 5; m++) {
	  forcing[i][j][k][m] = 0.0;
	}
      }
    }
  }
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
	buf[i][0] = cuf[i] + buf[i][2] * buf[i][2] + 
	  buf[i][3] * buf[i][3];
	q[i] = 0.5*(buf[i][1]*ue[i][1] + buf[i][2]*ue[i][2] +
		    buf[i][3]*ue[i][3]);
      }
      for (i = 1; i < grid_points[0]-1; i++) {
	im1 = i-1;
	ip1 = i+1;
	forcing[i][j][k][0] = forcing[i][j][k][0] -
	  tx2*(ue[ip1][1]-ue[im1][1])+
	  dx1tx1*(ue[ip1][0]-2.0*ue[i][0]+ue[im1][0]);
	forcing[i][j][k][1] = forcing[i][j][k][1] -
	  tx2 * ((ue[ip1][1]*buf[ip1][1]+c2*(ue[ip1][4]-q[ip1]))-
		 (ue[im1][1]*buf[im1][1]+c2*(ue[im1][4]-q[im1])))+
	  xxcon1*(buf[ip1][1]-2.0*buf[i][1]+buf[im1][1])+
	  dx2tx1*( ue[ip1][1]-2.0* ue[i][1]+ ue[im1][1]);
	forcing[i][j][k][2] = forcing[i][j][k][2] -
	  tx2 * (ue[ip1][2]*buf[ip1][1]-ue[im1][2]*buf[im1][1])+
	  xxcon2*(buf[ip1][2]-2.0*buf[i][2]+buf[im1][2])+
	  dx3tx1*( ue[ip1][2]-2.0* ue[i][2]+ ue[im1][2]);
	forcing[i][j][k][3] = forcing[i][j][k][3] -
	  tx2*(ue[ip1][3]*buf[ip1][1]-ue[im1][3]*buf[im1][1])+
	  xxcon2*(buf[ip1][3]-2.0*buf[i][3]+buf[im1][3])+
	  dx4tx1*( ue[ip1][3]-2.0* ue[i][3]+ ue[im1][3]);
	forcing[i][j][k][4] = forcing[i][j][k][4] -
	  tx2*(buf[ip1][1]*(c1*ue[ip1][4]-c2*q[ip1])-
	       buf[im1][1]*(c1*ue[im1][4]-c2*q[im1]))+
	  0.5*xxcon3*(buf[ip1][0]-2.0*buf[i][0]+buf[im1][0])+
	  xxcon4*(cuf[ip1]-2.0*cuf[i]+cuf[im1])+
	  xxcon5*(buf[ip1][4]-2.0*buf[i][4]+buf[im1][4])+
	  dx5tx1*( ue[ip1][4]-2.0* ue[i][4]+ ue[im1][4]);
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
	for (i = 1*3; i <= grid_points[0]-3*1-1; i++) {
	  forcing[i][j][k][m] = forcing[i][j][k][m] - dssp*
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
  for (i = 1; i < grid_points[0]-1; i++) {
    xi = (double)i * dnxm1;
    for (k = 1; k < grid_points[2]-1; k++) {
      zeta = (double)k * dnzm1;
      for (j = 0; j < grid_points[1]; j++) {
	eta = (double)j * dnym1;
	exact_solution(xi, eta, zeta, dtemp);
	for (m = 0; m < 5; m++) {
	  ue[j][m] = dtemp[m];
	}
	dtpp = 1.0/dtemp[0];
	for (m = 1; m <= 4; m++) {
	  buf[j][m] = dtpp * dtemp[m];
	}
	cuf[j]   = buf[j][2] * buf[j][2];
	buf[j][0] = cuf[j] + buf[j][1] * buf[j][1] + 
	  buf[j][3] * buf[j][3];
	q[j] = 0.5*(buf[j][1]*ue[j][1] + buf[j][2]*ue[j][2] +
		    buf[j][3]*ue[j][3]);
      }
      for (j = 1; j < grid_points[1]-1; j++) {
	jm1 = j-1;
	jp1 = j+1;
	forcing[i][j][k][0] = forcing[i][j][k][0] -
	  ty2*( ue[jp1][2]-ue[jm1][2] )+
	  dy1ty1*(ue[jp1][0]-2.0*ue[j][0]+ue[jm1][0]);
	forcing[i][j][k][1] = forcing[i][j][k][1] -
	  ty2*(ue[jp1][1]*buf[jp1][2]-ue[jm1][1]*buf[jm1][2])+
	  yycon2*(buf[jp1][1]-2.0*buf[j][1]+buf[jm1][1])+
	  dy2ty1*( ue[jp1][1]-2.0* ue[j][1]+ ue[jm1][1]);
	forcing[i][j][k][2] = forcing[i][j][k][2] -
	  ty2*((ue[jp1][2]*buf[jp1][2]+c2*(ue[jp1][4]-q[jp1]))-
	       (ue[jm1][2]*buf[jm1][2]+c2*(ue[jm1][4]-q[jm1])))+
	  yycon1*(buf[jp1][2]-2.0*buf[j][2]+buf[jm1][2])+
	  dy3ty1*( ue[jp1][2]-2.0*ue[j][2] +ue[jm1][2]);
	forcing[i][j][k][3] = forcing[i][j][k][3] -
	  ty2*(ue[jp1][3]*buf[jp1][2]-ue[jm1][3]*buf[jm1][2])+
	  yycon2*(buf[jp1][3]-2.0*buf[j][3]+buf[jm1][3])+
	  dy4ty1*( ue[jp1][3]-2.0*ue[j][3]+ ue[jm1][3]);
	forcing[i][j][k][4] = forcing[i][j][k][4] -
	  ty2*(buf[jp1][2]*(c1*ue[jp1][4]-c2*q[jp1])-
	       buf[jm1][2]*(c1*ue[jm1][4]-c2*q[jm1]))+
	  0.5*yycon3*(buf[jp1][0]-2.0*buf[j][0]+
                      buf[jm1][0])+
	  yycon4*(cuf[jp1]-2.0*cuf[j]+cuf[jm1])+
	  yycon5*(buf[jp1][4]-2.0*buf[j][4]+buf[jm1][4])+
	  dy5ty1*(ue[jp1][4]-2.0*ue[j][4]+ue[jm1][4]);
      }
      for (m = 0; m < 5; m++) {
	j = 1;
	forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
	  (5.0*ue[j][m] - 4.0*ue[j+1][m] +ue[j+2][m]);
	j = 2;
	forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
	  (-4.0*ue[j-1][m] + 6.0*ue[j][m] -
	   4.0*ue[j+1][m] +       ue[j+2][m]);
      }
      for (m = 0; m < 5; m++) {
	for (j = 1*3; j <= grid_points[1]-3*1-1; j++) {
	  forcing[i][j][k][m] = forcing[i][j][k][m] - dssp*
	    (ue[j-2][m] - 4.0*ue[j-1][m] +
	     6.0*ue[j][m] - 4.0*ue[j+1][m] + ue[j+2][m]);
	}
      }
      for (m = 0; m < 5; m++) {
	j = grid_points[1]-3;
	forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
	  (ue[j-2][m] - 4.0*ue[j-1][m] +
	   6.0*ue[j][m] - 4.0*ue[j+1][m]);
	j = grid_points[1]-2;
	forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
	  (ue[j-2][m] - 4.0*ue[j-1][m] + 5.0*ue[j][m]);
      }
    }
  }
  for (i = 1; i < grid_points[0]-1; i++) {
    xi = (double)i * dnxm1;
    for (j = 1; j < grid_points[1]-1; j++) {
      eta = (double)j * dnym1;
      for (k = 0; k < grid_points[2]; k++) {
	zeta = (double)k * dnzm1;
	exact_solution(xi, eta, zeta, dtemp);
	for (m = 0; m < 5; m++) {
	  ue[k][m] = dtemp[m];
	}
	dtpp = 1.0/dtemp[0];
	for (m = 1; m <= 4; m++) {
	  buf[k][m] = dtpp * dtemp[m];
	}
	cuf[k]   = buf[k][3] * buf[k][3];
	buf[k][0] = cuf[k] + buf[k][1] * buf[k][1] + 
	  buf[k][2] * buf[k][2];
	q[k] = 0.5*(buf[k][1]*ue[k][1] + buf[k][2]*ue[k][2] +
		    buf[k][3]*ue[k][3]);
      }
      for (k = 1; k < grid_points[2]-1; k++) {
	km1 = k-1;
	kp1 = k+1;
	forcing[i][j][k][0] = forcing[i][j][k][0] -
	  tz2*( ue[kp1][3]-ue[km1][3] )+
	  dz1tz1*(ue[kp1][0]-2.0*ue[k][0]+ue[km1][0]);
	forcing[i][j][k][1] = forcing[i][j][k][1] -
	  tz2 * (ue[kp1][1]*buf[kp1][3]-ue[km1][1]*buf[km1][3])+
	  zzcon2*(buf[kp1][1]-2.0*buf[k][1]+buf[km1][1])+
	  dz2tz1*( ue[kp1][1]-2.0* ue[k][1]+ ue[km1][1]);
	forcing[i][j][k][2] = forcing[i][j][k][2] -
	  tz2 * (ue[kp1][2]*buf[kp1][3]-ue[km1][2]*buf[km1][3])+
	  zzcon2*(buf[kp1][2]-2.0*buf[k][2]+buf[km1][2])+
	  dz3tz1*(ue[kp1][2]-2.0*ue[k][2]+ue[km1][2]);
	forcing[i][j][k][3] = forcing[i][j][k][3] -
	  tz2 * ((ue[kp1][3]*buf[kp1][3]+c2*(ue[kp1][4]-q[kp1]))-
		 (ue[km1][3]*buf[km1][3]+c2*(ue[km1][4]-q[km1])))+
	  zzcon1*(buf[kp1][3]-2.0*buf[k][3]+buf[km1][3])+
	  dz4tz1*( ue[kp1][3]-2.0*ue[k][3] +ue[km1][3]);
	forcing[i][j][k][4] = forcing[i][j][k][4] -
	  tz2 * (buf[kp1][3]*(c1*ue[kp1][4]-c2*q[kp1])-
		 buf[km1][3]*(c1*ue[km1][4]-c2*q[km1]))+
	  0.5*zzcon3*(buf[kp1][0]-2.0*buf[k][0]
                      +buf[km1][0])+
	  zzcon4*(cuf[kp1]-2.0*cuf[k]+cuf[km1])+
	  zzcon5*(buf[kp1][4]-2.0*buf[k][4]+buf[km1][4])+
	  dz5tz1*( ue[kp1][4]-2.0*ue[k][4]+ ue[km1][4]);
      }
      for (m = 0; m < 5; m++) {
	k = 1;
	forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
	  (5.0*ue[k][m] - 4.0*ue[k+1][m] +ue[k+2][m]);
	k = 2;
	forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
	  (-4.0*ue[k-1][m] + 6.0*ue[k][m] -
	   4.0*ue[k+1][m] +       ue[k+2][m]);
      }
      for (m = 0; m < 5; m++) {
	for (k = 1*3; k <= grid_points[2]-3*1-1; k++) {
	  forcing[i][j][k][m] = forcing[i][j][k][m] - dssp*
	    (ue[k-2][m] - 4.0*ue[k-1][m] +
	     6.0*ue[k][m] - 4.0*ue[k+1][m] + ue[k+2][m]);
	}
      }
      for (m = 0; m < 5; m++) {
	k = grid_points[2]-3;
	forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
	  (ue[k-2][m] - 4.0*ue[k-1][m] +
	   6.0*ue[k][m] - 4.0*ue[k+1][m]);
	k = grid_points[2]-2;
	forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
	  (ue[k-2][m] - 4.0*ue[k-1][m] + 5.0*ue[k][m]);
      }
    }
  }
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	for (m = 0; m < 5; m++) {
	  forcing[i][j][k][m] = -1.0 * forcing[i][j][k][m];
	}
      }
    }
  }
}
}
static void exact_solution(double xi, double eta, double zeta,
			   double dtemp[5]) {
  int m;
  for (m = 0; m < 5; m++) {
    dtemp[m] =  ce[m][0] +
      xi*(ce[m][1] + xi*(ce[m][4] + xi*(ce[m][7]
					+ xi*ce[m][10]))) +
      eta*(ce[m][2] + eta*(ce[m][5] + eta*(ce[m][8]
					   + eta*ce[m][11])))+
      zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + 
					      zeta*ce[m][12])));
  }
}
static void initialize(void) {
{
  int i, j, k, m, ix, iy, iz;
  double xi, eta, zeta, Pface[2][3][5], Pxi, Peta, Pzeta, temp[5];
  for (i = 0; i < IMAX; i++) {
    for (j = 0; j < IMAX; j++) {
      for (k = 0; k < IMAX; k++) {
	for (m = 0; m < 5; m++) {
	  u[i][j][k][m] = 1.0;
	}
      }
    }
  }
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      for (k = 0; k < grid_points[2]; k++) {
	zeta = (double)k * dnzm1;
	for (ix = 0; ix < 2; ix++) {
	  exact_solution((double)ix, eta, zeta, 
                         &(Pface[ix][0][0]));
	}
	for (iy = 0; iy < 2; iy++) {
	  exact_solution(xi, (double)iy , zeta, 
                         &Pface[iy][1][0]);
	}
	for (iz = 0; iz < 2; iz++) {
	  exact_solution(xi, eta, (double)iz,   
                         &Pface[iz][2][0]);
	}
	for (m = 0; m < 5; m++) {
	  Pxi   = xi   * Pface[1][0][m] + 
	    (1.0-xi)   * Pface[0][0][m];
	  Peta  = eta  * Pface[1][1][m] + 
	    (1.0-eta)  * Pface[0][1][m];
	  Pzeta = zeta * Pface[1][2][m] + 
	    (1.0-zeta) * Pface[0][2][m];
	  u[i][j][k][m] = Pxi + Peta + Pzeta - 
	    Pxi*Peta - Pxi*Pzeta - Peta*Pzeta + 
	    Pxi*Peta*Pzeta;
	}
      }
    }
  }
  i = 0;
  xi = 0.0;
  for (j = 0; j < grid_points[1]; j++) {
    eta = (double)j * dnym1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m];
      }
    }
  }
  i = grid_points[0]-1;
  xi = 1.0;
  for (j = 0; j < grid_points[1]; j++) {
    eta = (double)j * dnym1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m];
      }
    }
  }
  j = 0;
  eta = 0.0;
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m];
      }
    }
  }
  j = grid_points[1]-1;
  eta = 1.0;
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (k = 0; k < grid_points[2]; k++) {
      zeta = (double)k * dnzm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m];
      }
    }
  }
  k = 0;
  zeta = 0.0;
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i *dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m];
      }
    }
  }
  k = grid_points[2]-1;
  zeta = 1.0;
  for (i = 0; i < grid_points[0]; i++) {
    xi = (double)i * dnxm1;
    for (j = 0; j < grid_points[1]; j++) {
      eta = (double)j * dnym1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
	u[i][j][k][m] = temp[m];
      }
    }
  }
}
}
static void lhsinit(void) {
{
  int i, j, k, m, n;
  for (i = 0; i < grid_points[0]; i++) {
    for (j = 0; j < grid_points[1]; j++) {
      for (k = 0; k < grid_points[2]; k++) {
	for (m = 0; m < 5; m++) {
	  for (n = 0; n < 5; n++) {
	    lhs[i][j][k][0][m][n] = 0.0;
	    lhs[i][j][k][1][m][n] = 0.0;
	    lhs[i][j][k][2][m][n] = 0.0;
	  }
	}
      }
    }
  }
  for (i = 0; i < grid_points[0]; i++) {
    for (j = 0; j < grid_points[1]; j++) {
      for (k = 0; k < grid_points[2]; k++) {
	for (m = 0; m < 5; m++) {
	  lhs[i][j][k][1][m][m] = 1.0;
	}
      }
    }
  }
}
}
#include <omp.h> // 使用 OpenMP 魔法，需要先告诉电脑

static void lhsx(void) {
  int i, j, k;
  // 这几个是临时变量，就像是临时用的小碗
  // 把它们声明在循环里面，或者设为 private，确保每个厨师有自己的碗
  double tmp1, tmp2, tmp3; 

  // 魔法开始！把所有小厨师叫过来
  // private(...)：括号里的是每个小厨师自己独有的工具
  #pragma omp parallel private(i, j, k, tmp1, tmp2, tmp3)
  {
    // #pragma omp for：把下面的任务分给不同的小厨师
    // collapse(2)：把 j 和 k 两层循环的所有组合任务，打包分给大家
    #pragma omp for collapse(2)
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
        // 对于一个固定的 (j,k) 任务，下面的步骤是连贯的，由一个厨师独立完成

        // 第1步: 准备基础材料 (fjac 和 njac)
        for (i = 0; i < grid_points[0]; i++) {
          tmp1 = 1.0 / u[i][j][k][0];
          tmp2 = tmp1 * tmp1;
          tmp3 = tmp1 * tmp2;
          // ... (省略了中间一长串一模一样的 fjac, njac 计算)
          fjac[ i][ j][ k][0][0] = 0.0;
          fjac[ i][ j][ k][0][1] = 1.0;
          fjac[ i][ j][ k][0][2] = 0.0;
          fjac[ i][ j][ k][0][3] = 0.0;
          fjac[ i][ j][ k][0][4] = 0.0;
          fjac[ i][ j][ k][1][0] = -(u[i][j][k][1] * tmp2 * u[i][j][k][1]) + c2 * 0.50 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3] ) * tmp2;
          fjac[i][j][k][1][1] = ( 2.0 - c2 ) * ( u[i][j][k][1] / u[i][j][k][0] );
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
          fjac[i][j][k][4][0] = ( c2 * ( u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3] ) * tmp2 - c1 * ( u[i][j][k][4] * tmp1 ) ) * ( u[i][j][k][1] * tmp1 );
          fjac[i][j][k][4][1] = c1 * u[i][j][k][4] * tmp1 - 0.50 * c2 * ( 3.0*u[i][j][k][1]*u[i][j][k][1] + u[i][j][k][2]*u[i][j][k][2] + u[i][j][k][3]*u[i][j][k][3] ) * tmp2;
          fjac[i][j][k][4][2] = - c2 * ( u[i][j][k][2]*u[i][j][k][1] ) * tmp2;
          fjac[i][j][k][4][3] = - c2 * ( u[i][j][k][3]*u[i][j][k][1] ) * tmp2;
          fjac[i][j][k][4][4] = c1 * ( u[i][j][k][1] * tmp1 );
          njac[i][j][k][0][0] = 0.0;
          njac[i][j][k][0][1] = 0.0;
          njac[i][j][k][0][2] = 0.0;
          njac[i][j][k][0][3] = 0.0;
          njac[i][j][k][0][4] = 0.0;
          njac[i][j][k][1][0] = - con43 * c3c4 * tmp2 * u[i][j][k][1];
          njac[i][j][k][1][1] = con43 * c3c4 * tmp1;
          njac[i][j][k][1][2] = 0.0;
          njac[i][j][k][1][3] = 0.0;
          njac[i][j][k][1][4] = 0.0;
          njac[i][j][k][2][0] = - c3c4 * tmp2 * u[i][j][k][2];
          njac[i][j][k][2][1] = 0.0;
          njac[i][j][k][2][2] = c3c4 * tmp1;
          njac[i][j][k][2][3] = 0.0;
          njac[i][j][k][2][4] = 0.0;
          njac[i][j][k][3][0] = - c3c4 * tmp2 * u[i][j][k][3];
          njac[i][j][k][3][1] = 0.0;
          njac[i][j][k][3][2] = 0.0;
          njac[i][j][k][3][3] = c3c4 * tmp1;
          njac[i][j][k][3][4] = 0.0;
          njac[i][j][k][4][0] = - ( con43 * c3c4 - c1345 ) * tmp3 * (u[i][j][k][1]*u[i][j][k][1]) - ( c3c4 - c1345 ) * tmp3 * (u[i][j][k][2]*u[i][j][k][2]) - ( c3c4 - c1345 ) * tmp3 * (u[i][j][k][3]*u[i][j][k][3]) - c1345 * tmp2 * u[i][j][k][4];
          njac[i][j][k][4][1] = ( con43 * c3c4 - c1345 ) * tmp2 * u[i][j][k][1];
          njac[i][j][k][4][2] = ( c3c4 - c1345 ) * tmp2 * u[i][j][k][2];
          njac[i][j][k][4][3] = ( c3c4 - c1345 ) * tmp2 * u[i][j][k][3];
          njac[i][j][k][4][4] = ( c1345 ) * tmp1;
        }

        // 第2步: 用基础材料调配酱料 (lhs)
        for (i = 1; i < grid_points[0]-1; i++) {
          tmp1 = dt * tx1;
          tmp2 = dt * tx2;
          // ... (省略了中间一长串一模一样的 lhs 计算)
          lhs[i][j][k][AA][0][0] = - tmp2 * fjac[i-1][j][k][0][0] - tmp1 * njac[i-1][j][k][0][0] - tmp1 * dx1;
          lhs[i][j][k][AA][0][1] = - tmp2 * fjac[i-1][j][k][0][1] - tmp1 * njac[i-1][j][k][0][1];
          lhs[i][j][k][AA][0][2] = - tmp2 * fjac[i-1][j][k][0][2] - tmp1 * njac[i-1][j][k][0][2];
          lhs[i][j][k][AA][0][3] = - tmp2 * fjac[i-1][j][k][0][3] - tmp1 * njac[i-1][j][k][0][3];
          lhs[i][j][k][AA][0][4] = - tmp2 * fjac[i-1][j][k][0][4] - tmp1 * njac[i-1][j][k][0][4];
          lhs[i][j][k][AA][1][0] = - tmp2 * fjac[i-1][j][k][1][0] - tmp1 * njac[i-1][j][k][1][0];
          lhs[i][j][k][AA][1][1] = - tmp2 * fjac[i-1][j][k][1][1] - tmp1 * njac[i-1][j][k][1][1] - tmp1 * dx2;
          lhs[i][j][k][AA][1][2] = - tmp2 * fjac[i-1][j][k][1][2] - tmp1 * njac[i-1][j][k][1][2];
          lhs[i][j][k][AA][1][3] = - tmp2 * fjac[i-1][j][k][1][3] - tmp1 * njac[i-1][j][k][1][3];
          lhs[i][j][k][AA][1][4] = - tmp2 * fjac[i-1][j][k][1][4] - tmp1 * njac[i-1][j][k][1][4];
          lhs[i][j][k][AA][2][0] = - tmp2 * fjac[i-1][j][k][2][0] - tmp1 * njac[i-1][j][k][2][0];
          lhs[i][j][k][AA][2][1] = - tmp2 * fjac[i-1][j][k][2][1] - tmp1 * njac[i-1][j][k][2][1];
          lhs[i][j][k][AA][2][2] = - tmp2 * fjac[i-1][j][k][2][2] - tmp1 * njac[i-1][j][k][2][2] - tmp1 * dx3;
          lhs[i][j][k][AA][2][3] = - tmp2 * fjac[i-1][j][k][2][3] - tmp1 * njac[i-1][j][k][2][3];
          lhs[i][j][k][AA][2][4] = - tmp2 * fjac[i-1][j][k][2][4] - tmp1 * njac[i-1][j][k][2][4];
          lhs[i][j][k][AA][3][0] = - tmp2 * fjac[i-1][j][k][3][0] - tmp1 * njac[i-1][j][k][3][0];
          lhs[i][j][k][AA][3][1] = - tmp2 * fjac[i-1][j][k][3][1] - tmp1 * njac[i-1][j][k][3][1];
          lhs[i][j][k][AA][3][2] = - tmp2 * fjac[i-1][j][k][3][2] - tmp1 * njac[i-1][j][k][3][2];
          lhs[i][j][k][AA][3][3] = - tmp2 * fjac[i-1][j][k][3][3] - tmp1 * njac[i-1][j][k][3][3] - tmp1 * dx4;
          lhs[i][j][k][AA][3][4] = - tmp2 * fjac[i-1][j][k][3][4] - tmp1 * njac[i-1][j][k][3][4];
          lhs[i][j][k][AA][4][0] = - tmp2 * fjac[i-1][j][k][4][0] - tmp1 * njac[i-1][j][k][4][0];
          lhs[i][j][k][AA][4][1] = - tmp2 * fjac[i-1][j][k][4][1] - tmp1 * njac[i-1][j][k][4][1];
          lhs[i][j][k][AA][4][2] = - tmp2 * fjac[i-1][j][k][4][2] - tmp1 * njac[i-1][j][k][4][2];
          lhs[i][j][k][AA][4][3] = - tmp2 * fjac[i-1][j][k][4][3] - tmp1 * njac[i-1][j][k][4][3];
          lhs[i][j][k][AA][4][4] = - tmp2 * fjac[i-1][j][k][4][4] - tmp1 * njac[i-1][j][k][4][4] - tmp1 * dx5;
          lhs[i][j][k][BB][0][0] = 1.0 + tmp1 * 2.0 * njac[i][j][k][0][0] + tmp1 * 2.0 * dx1;
          lhs[i][j][k][BB][0][1] = tmp1 * 2.0 * njac[i][j][k][0][1];
          lhs[i][j][k][BB][0][2] = tmp1 * 2.0 * njac[i][j][k][0][2];
          lhs[i][j][k][BB][0][3] = tmp1 * 2.0 * njac[i][j][k][0][3];
          lhs[i][j][k][BB][0][4] = tmp1 * 2.0 * njac[i][j][k][0][4];
          lhs[i][j][k][BB][1][0] = tmp1 * 2.0 * njac[i][j][k][1][0];
          lhs[i][j][k][BB][1][1] = 1.0 + tmp1 * 2.0 * njac[i][j][k][1][1] + tmp1 * 2.0 * dx2;
          lhs[i][j][k][BB][1][2] = tmp1 * 2.0 * njac[i][j][k][1][2];
          lhs[i][j][k][BB][1][3] = tmp1 * 2.0 * njac[i][j][k][1][3];
          lhs[i][j][k][BB][1][4] = tmp1 * 2.0 * njac[i][j][k][1][4];
          lhs[i][j][k][BB][2][0] = tmp1 * 2.0 * njac[i][j][k][2][0];
          lhs[i][j][k][BB][2][1] = tmp1 * 2.0 * njac[i][j][k][2][1];
          lhs[i][j][k][BB][2][2] = 1.0 + tmp1 * 2.0 * njac[i][j][k][2][2] + tmp1 * 2.0 * dx3;
          lhs[i][j][k][BB][2][3] = tmp1 * 2.0 * njac[i][j][k][2][3];
          lhs[i][j][k][BB][2][4] = tmp1 * 2.0 * njac[i][j][k][2][4];
          lhs[i][j][k][BB][3][0] = tmp1 * 2.0 * njac[i][j][k][3][0];
          lhs[i][j][k][BB][3][1] = tmp1 * 2.0 * njac[i][j][k][3][1];
          lhs[i][j][k][BB][3][2] = tmp1 * 2.0 * njac[i][j][k][3][2];
          lhs[i][j][k][BB][3][3] = 1.0 + tmp1 * 2.0 * njac[i][j][k][3][3] + tmp1 * 2.0 * dx4;
          lhs[i][j][k][BB][3][4] = tmp1 * 2.0 * njac[i][j][k][3][4];
          lhs[i][j][k][BB][4][0] = tmp1 * 2.0 * njac[i][j][k][4][0];
          lhs[i][j][k][BB][4][1] = tmp1 * 2.0 * njac[i][j][k][4][1];
          lhs[i][j][k][BB][4][2] = tmp1 * 2.0 * njac[i][j][k][4][2];
          lhs[i][j][k][BB][4][3] = tmp1 * 2.0 * njac[i][j][k][4][3];
          lhs[i][j][k][BB][4][4] = 1.0 + tmp1 * 2.0 * njac[i][j][k][4][4] + tmp1 * 2.0 * dx5;
          lhs[i][j][k][CC][0][0] = tmp2 * fjac[i+1][j][k][0][0] - tmp1 * njac[i+1][j][k][0][0] - tmp1 * dx1;
          lhs[i][j][k][CC][0][1] = tmp2 * fjac[i+1][j][k][0][1] - tmp1 * njac[i+1][j][k][0][1];
          lhs[i][j][k][CC][0][2] = tmp2 * fjac[i+1][j][k][0][2] - tmp1 * njac[i+1][j][k][0][2];
          lhs[i][j][k][CC][0][3] = tmp2 * fjac[i+1][j][k][0][3] - tmp1 * njac[i+1][j][k][0][3];
          lhs[i][j][k][CC][0][4] = tmp2 * fjac[i+1][j][k][0][4] - tmp1 * njac[i+1][j][k][0][4];
          lhs[i][j][k][CC][1][0] = tmp2 * fjac[i+1][j][k][1][0] - tmp1 * njac[i+1][j][k][1][0];
          lhs[i][j][k][CC][1][1] = tmp2 * fjac[i+1][j][k][1][1] - tmp1 * njac[i+1][j][k][1][1] - tmp1 * dx2;
          lhs[i][j][k][CC][1][2] = tmp2 * fjac[i+1][j][k][1][2] - tmp1 * njac[i+1][j][k][1][2];
          lhs[i][j][k][CC][1][3] = tmp2 * fjac[i+1][j][k][1][3] - tmp1 * njac[i+1][j][k][1][3];
          lhs[i][j][k][CC][1][4] = tmp2 * fjac[i+1][j][k][1][4] - tmp1 * njac[i+1][j][k][1][4];
          lhs[i][j][k][CC][2][0] = tmp2 * fjac[i+1][j][k][2][0] - tmp1 * njac[i+1][j][k][2][0];
          lhs[i][j][k][CC][2][1] = tmp2 * fjac[i+1][j][k][2][1] - tmp1 * njac[i+1][j][k][2][1];
          lhs[i][j][k][CC][2][2] = tmp2 * fjac[i+1][j][k][2][2] - tmp1 * njac[i+1][j][k][2][2] - tmp1 * dx3;
          lhs[i][j][k][CC][2][3] = tmp2 * fjac[i+1][j][k][2][3] - tmp1 * njac[i+1][j][k][2][3];
          lhs[i][j][k][CC][2][4] = tmp2 * fjac[i+1][j][k][2][4] - tmp1 * njac[i+1][j][k][2][4];
          lhs[i][j][k][CC][3][0] = tmp2 * fjac[i+1][j][k][3][0] - tmp1 * njac[i+1][j][k][3][0];
          lhs[i][j][k][CC][3][1] = tmp2 * fjac[i+1][j][k][3][1] - tmp1 * njac[i+1][j][k][3][1];
          lhs[i][j][k][CC][3][2] = tmp2 * fjac[i+1][j][k][3][2] - tmp1 * njac[i+1][j][k][3][2];
          lhs[i][j][k][CC][3][3] = tmp2 * fjac[i+1][j][k][3][3] - tmp1 * njac[i+1][j][k][3][3] - tmp1 * dx4;
          lhs[i][j][k][CC][3][4] = tmp2 * fjac[i+1][j][k][3][4] - tmp1 * njac[i+1][j][k][3][4];
          lhs[i][j][k][CC][4][0] = tmp2 * fjac[i+1][j][k][4][0] - tmp1 * njac[i+1][j][k][4][0];
          lhs[i][j][k][CC][4][1] = tmp2 * fjac[i+1][j][k][4][1] - tmp1 * njac[i+1][j][k][4][1];
          lhs[i][j][k][CC][4][2] = tmp2 * fjac[i+1][j][k][4][2] - tmp1 * njac[i+1][j][k][4][2];
          lhs[i][j][k][CC][4][3] = tmp2 * fjac[i+1][j][k][4][3] - tmp1 * njac[i+1][j][k][4][3];
          lhs[i][j][k][CC][4][4] = tmp2 * fjac[i+1][j][k][4][4] - tmp1 * njac[i+1][j][k][4][4] - tmp1 * dx5;
        }
      }
    }
  } // 团队合作区结束
}
#include <omp.h> // 引入OpenMP头文件
#include <math.h> // 为了使用 pow2，假设它是 pow(x, 2) 的宏或函数

// 假设这些变量和宏在其他地方已经定义
// extern double u[...][...][...][...];
// extern double fjac[...][...][...][...][...];
// extern double njac[...][...][...][...][...];
// extern double lhs[...][...][...][...][...][...];
// extern int grid_points[3];
// extern double c1, c2, c3c4, con43, c1345;
// extern double dt, ty1, ty2, dy1, dy2, dy3, dy4, dy5;
// #define AA 0 // 假设的宏定义
// #define BB 1
// #define CC 2
// #define pow2(x) ((x)*(x)) // 假设 pow2 是这样实现的

static void lhsy(void) {
  int i, j, k;
  double tmp1, tmp2, tmp3;

  // 仅创建一个并行区域，以满足约束
  // 将循环中使用的临时变量声明为线程私有，避免数据竞争
  #pragma omp parallel private(i, j, k, tmp1, tmp2, tmp3)
  {
    // =====================================================================
    // 第一个循环：计算 fjac 和 njac 矩阵
    // 使用 collapse(2) 将外两层循环合并，以改善负载均衡
    // 使用静态调度，因为每次迭代的工作量相似
    // =====================================================================
    #pragma omp for collapse(2) schedule(static)
    for (i = 1; i < grid_points[0]-1; i++) {
      for (j = 0; j < grid_points[1]; j++) {
        for (k = 1; k < grid_points[2]-1; k++) {
          tmp1 = 1.0 / u[i][j][k][0];
          tmp2 = tmp1 * tmp1;
          tmp3 = tmp1 * tmp2;

          fjac[ i][ j][ k][0][0] = 0.0;
          fjac[ i][ j][ k][0][1] = 0.0;
          fjac[ i][ j][ k][0][2] = 1.0;
          fjac[ i][ j][ k][0][3] = 0.0;
          fjac[ i][ j][ k][0][4] = 0.0;

          fjac[i][j][k][1][0] = - ( u[i][j][k][1]*u[i][j][k][2] )
            * tmp2;
          fjac[i][j][k][1][1] = u[i][j][k][2] * tmp1;
          fjac[i][j][k][1][2] = u[i][j][k][1] * tmp1;
          fjac[i][j][k][1][3] = 0.0;
          fjac[i][j][k][1][4] = 0.0;

          fjac[i][j][k][2][0] = - ( u[i][j][k][2]*u[i][j][k][2]*tmp2)
            + 0.50 * c2 * ( (  u[i][j][k][1] * u[i][j][k][1]
                               + u[i][j][k][2] * u[i][j][k][2]
                               + u[i][j][k][3] * u[i][j][k][3] )
                            * tmp2 );
          fjac[i][j][k][2][1] = - c2 *  u[i][j][k][1] * tmp1;
          fjac[i][j][k][2][2] = ( 2.0 - c2 )
            *  u[i][j][k][2] * tmp1;
          fjac[i][j][k][2][3] = - c2 * u[i][j][k][3] * tmp1;
          fjac[i][j][k][2][4] = c2;

          fjac[i][j][k][3][0] = - ( u[i][j][k][2]*u[i][j][k][3] )
            * tmp2;
          fjac[i][j][k][3][1] = 0.0;
          fjac[i][j][k][3][2] = u[i][j][k][3] * tmp1;
          fjac[i][j][k][3][3] = u[i][j][k][2] * tmp1;
          fjac[i][j][k][3][4] = 0.0;

          fjac[i][j][k][4][0] = ( c2 * (  u[i][j][k][1] * u[i][j][k][1]
                                          + u[i][j][k][2] * u[i][j][k][2]
                                          + u[i][j][k][3] * u[i][j][k][3] )
                                  * tmp2
                                  - c1 * u[i][j][k][4] * tmp1 ) 
            * u[i][j][k][2] * tmp1;
          fjac[i][j][k][4][1] = - c2 * u[i][j][k][1]*u[i][j][k][2] 
            * tmp2;
          fjac[i][j][k][4][2] = c1 * u[i][j][k][4] * tmp1 
            - 0.50 * c2 
            * ( (  u[i][j][k][1]*u[i][j][k][1]
                   + 3.0 * u[i][j][k][2]*u[i][j][k][2]
                   + u[i][j][k][3]*u[i][j][k][3] )
                * tmp2 );
          fjac[i][j][k][4][3] = - c2 * ( u[i][j][k][2]*u[i][j][k][3] )
            * tmp2;
          fjac[i][j][k][4][4] = c1 * u[i][j][k][2] * tmp1; 

          njac[i][j][k][0][0] = 0.0;
          njac[i][j][k][0][1] = 0.0;
          njac[i][j][k][0][2] = 0.0;
          njac[i][j][k][0][3] = 0.0;
          njac[i][j][k][0][4] = 0.0;

          njac[i][j][k][1][0] = - c3c4 * tmp2 * u[i][j][k][1];
          njac[i][j][k][1][1] =   c3c4 * tmp1;
          njac[i][j][k][1][2] =   0.0;
          njac[i][j][k][1][3] =   0.0;
          njac[i][j][k][1][4] =   0.0;

          njac[i][j][k][2][0] = - con43 * c3c4 * tmp2 * u[i][j][k][2];
          njac[i][j][k][2][1] =   0.0;
          njac[i][j][k][2][2] =   con43 * c3c4 * tmp1;
          njac[i][j][k][2][3] =   0.0;
          njac[i][j][k][2][4] =   0.0;

          njac[i][j][k][3][0] = - c3c4 * tmp2 * u[i][j][k][3];
          njac[i][j][k][3][1] =   0.0;
          njac[i][j][k][3][2] =   0.0;
          njac[i][j][k][3][3] =   c3c4 * tmp1;
          njac[i][j][k][3][4] =   0.0;

          njac[i][j][k][4][0] = - (  c3c4
                - c1345 ) * tmp3 * (pow2(u[i][j][k][1]))
            - ( con43 * c3c4
                - c1345 ) * tmp3 * (pow2(u[i][j][k][2]))
            - ( c3c4 - c1345 ) * tmp3 * (pow2(u[i][j][k][3]))
            - c1345 * tmp2 * u[i][j][k][4];
          njac[i][j][k][4][1] = (  c3c4 - c1345 ) * tmp2 * u[i][j][k][1];
          njac[i][j][k][4][2] = ( con43 * c3c4
                                  - c1345 ) * tmp2 * u[i][j][k][2];
          njac[i][j][k][4][3] = ( c3c4 - c1345 ) * tmp2 * u[i][j][k][3];
          njac[i][j][k][4][4] = ( c1345 ) * tmp1;
        }
      }
    }

    // 在此位置有一个隐式的屏障(barrier)。
    // 所有线程必须完成上面的 #pragma omp for 循环后，才能继续执行下面的代码。
    // 这保证了在计算 lhs 之前，所有的 fjac 和 njac 都已计算完毕。

    // =====================================================================
    // 第二个循环：计算 lhs 矩阵
    // 同样，这个循环的迭代之间没有依赖，可以安全地并行化。
    // =====================================================================
    #pragma omp for collapse(2) schedule(static)
    for (i = 1; i < grid_points[0]-1; i++) {
      for (j = 1; j < grid_points[1]-1; j++) {
        for (k = 1; k < grid_points[2]-1; k++) {
          tmp1 = dt * ty1;
          tmp2 = dt * ty2;
          lhs[i][j][k][AA][0][0] = - tmp2 * fjac[i][j-1][k][0][0]
            - tmp1 * njac[i][j-1][k][0][0]
            - tmp1 * dy1;
          lhs[i][j][k][AA][0][1] = - tmp2 * fjac[i][j-1][k][0][1]
            - tmp1 * njac[i][j-1][k][0][1];
          lhs[i][j][k][AA][0][2] = - tmp2 * fjac[i][j-1][k][0][2]
            - tmp1 * njac[i][j-1][k][0][2];
          lhs[i][j][k][AA][0][3] = - tmp2 * fjac[i][j-1][k][0][3]
            - tmp1 * njac[i][j-1][k][0][3];
          lhs[i][j][k][AA][0][4] = - tmp2 * fjac[i][j-1][k][0][4]
            - tmp1 * njac[i][j-1][k][0][4];
          lhs[i][j][k][AA][1][0] = - tmp2 * fjac[i][j-1][k][1][0]
            - tmp1 * njac[i][j-1][k][1][0];
          lhs[i][j][k][AA][1][1] = - tmp2 * fjac[i][j-1][k][1][1]
            - tmp1 * njac[i][j-1][k][1][1]
            - tmp1 * dy2;
          lhs[i][j][k][AA][1][2] = - tmp2 * fjac[i][j-1][k][1][2]
            - tmp1 * njac[i][j-1][k][1][2];
          lhs[i][j][k][AA][1][3] = - tmp2 * fjac[i][j-1][k][1][3]
            - tmp1 * njac[i][j-1][k][1][3];
          lhs[i][j][k][AA][1][4] = - tmp2 * fjac[i][j-1][k][1][4]
            - tmp1 * njac[i][j-1][k][1][4];
          lhs[i][j][k][AA][2][0] = - tmp2 * fjac[i][j-1][k][2][0]
            - tmp1 * njac[i][j-1][k][2][0];
          lhs[i][j][k][AA][2][1] = - tmp2 * fjac[i][j-1][k][2][1]
            - tmp1 * njac[i][j-1][k][2][1];
          lhs[i][j][k][AA][2][2] = - tmp2 * fjac[i][j-1][k][2][2]
            - tmp1 * njac[i][j-1][k][2][2]
            - tmp1 * dy3;
          lhs[i][j][k][AA][2][3] = - tmp2 * fjac[i][j-1][k][2][3]
            - tmp1 * njac[i][j-1][k][2][3];
          lhs[i][j][k][AA][2][4] = - tmp2 * fjac[i][j-1][k][2][4]
            - tmp1 * njac[i][j-1][k][2][4];
          lhs[i][j][k][AA][3][0] = - tmp2 * fjac[i][j-1][k][3][0]
            - tmp1 * njac[i][j-1][k][3][0];
          lhs[i][j][k][AA][3][1] = - tmp2 * fjac[i][j-1][k][3][1]
            - tmp1 * njac[i][j-1][k][3][1];
          lhs[i][j][k][AA][3][2] = - tmp2 * fjac[i][j-1][k][3][2]
            - tmp1 * njac[i][j-1][k][3][2];
          lhs[i][j][k][AA][3][3] = - tmp2 * fjac[i][j-1][k][3][3]
            - tmp1 * njac[i][j-1][k][3][3]
            - tmp1 * dy4;
          lhs[i][j][k][AA][3][4] = - tmp2 * fjac[i][j-1][k][3][4]
            - tmp1 * njac[i][j-1][k][3][4];
          lhs[i][j][k][AA][4][0] = - tmp2 * fjac[i][j-1][k][4][0]
            - tmp1 * njac[i][j-1][k][4][0];
          lhs[i][j][k][AA][4][1] = - tmp2 * fjac[i][j-1][k][4][1]
            - tmp1 * njac[i][j-1][k][4][1];
          lhs[i][j][k][AA][4][2] = - tmp2 * fjac[i][j-1][k][4][2]
            - tmp1 * njac[i][j-1][k][4][2];
          lhs[i][j][k][AA][4][3] = - tmp2 * fjac[i][j-1][k][4][3]
            - tmp1 * njac[i][j-1][k][4][3];
          lhs[i][j][k][AA][4][4] = - tmp2 * fjac[i][j-1][k][4][4]
            - tmp1 * njac[i][j-1][k][4][4]
            - tmp1 * dy5;
          lhs[i][j][k][BB][0][0] = 1.0
            + tmp1 * 2.0 * njac[i][j][k][0][0]
            + tmp1 * 2.0 * dy1;
          lhs[i][j][k][BB][0][1] = tmp1 * 2.0 * njac[i][j][k][0][1];
          lhs[i][j][k][BB][0][2] = tmp1 * 2.0 * njac[i][j][k][0][2];
          lhs[i][j][k][BB][0][3] = tmp1 * 2.0 * njac[i][j][k][0][3];
          lhs[i][j][k][BB][0][4] = tmp1 * 2.0 * njac[i][j][k][0][4];
          lhs[i][j][k][BB][1][0] = tmp1 * 2.0 * njac[i][j][k][1][0];
          lhs[i][j][k][BB][1][1] = 1.0
            + tmp1 * 2.0 * njac[i][j][k][1][1]
            + tmp1 * 2.0 * dy2;
          lhs[i][j][k][BB][1][2] = tmp1 * 2.0 * njac[i][j][k][1][2];
          lhs[i][j][k][BB][1][3] = tmp1 * 2.0 * njac[i][j][k][1][3];
          lhs[i][j][k][BB][1][4] = tmp1 * 2.0 * njac[i][j][k][1][4];
          lhs[i][j][k][BB][2][0] = tmp1 * 2.0 * njac[i][j][k][2][0];
          lhs[i][j][k][BB][2][1] = tmp1 * 2.0 * njac[i][j][k][2][1];
          lhs[i][j][k][BB][2][2] = 1.0
            + tmp1 * 2.0 * njac[i][j][k][2][2]
            + tmp1 * 2.0 * dy3;
          lhs[i][j][k][BB][2][3] = tmp1 * 2.0 * njac[i][j][k][2][3];
          lhs[i][j][k][BB][2][4] = tmp1 * 2.0 * njac[i][j][k][2][4];
          lhs[i][j][k][BB][3][0] = tmp1 * 2.0 * njac[i][j][k][3][0];
          lhs[i][j][k][BB][3][1] = tmp1 * 2.0 * njac[i][j][k][3][1];
          lhs[i][j][k][BB][3][2] = tmp1 * 2.0 * njac[i][j][k][3][2];
          lhs[i][j][k][BB][3][3] = 1.0
            + tmp1 * 2.0 * njac[i][j][k][3][3]
            + tmp1 * 2.0 * dy4;
          lhs[i][j][k][BB][3][4] = tmp1 * 2.0 * njac[i][j][k][3][4];
          lhs[i][j][k][BB][4][0] = tmp1 * 2.0 * njac[i][j][k][4][0];
          lhs[i][j][k][BB][4][1] = tmp1 * 2.0 * njac[i][j][k][4][1];
          lhs[i][j][k][BB][4][2] = tmp1 * 2.0 * njac[i][j][k][4][2];
          lhs[i][j][k][BB][4][3] = tmp1 * 2.0 * njac[i][j][k][4][3];
          lhs[i][j][k][BB][4][4] = 1.0
            + tmp1 * 2.0 * njac[i][j][k][4][4] 
            + tmp1 * 2.0 * dy5;
          lhs[i][j][k][CC][0][0] =  tmp2 * fjac[i][j+1][k][0][0]
            - tmp1 * njac[i][j+1][k][0][0]
            - tmp1 * dy1;
          lhs[i][j][k][CC][0][1] =  tmp2 * fjac[i][j+1][k][0][1]
            - tmp1 * njac[i][j+1][k][0][1];
          lhs[i][j][k][CC][0][2] =  tmp2 * fjac[i][j+1][k][0][2]
            - tmp1 * njac[i][j+1][k][0][2];
          lhs[i][j][k][CC][0][3] =  tmp2 * fjac[i][j+1][k][0][3]
            - tmp1 * njac[i][j+1][k][0][3];
          lhs[i][j][k][CC][0][4] =  tmp2 * fjac[i][j+1][k][0][4]
            - tmp1 * njac[i][j+1][k][0][4];
          lhs[i][j][k][CC][1][0] =  tmp2 * fjac[i][j+1][k][1][0]
            - tmp1 * njac[i][j+1][k][1][0];
          lhs[i][j][k][CC][1][1] =  tmp2 * fjac[i][j+1][k][1][1]
            - tmp1 * njac[i][j+1][k][1][1]
            - tmp1 * dy2;
          lhs[i][j][k][CC][1][2] =  tmp2 * fjac[i][j+1][k][1][2]
            - tmp1 * njac[i][j+1][k][1][2];
          lhs[i][j][k][CC][1][3] =  tmp2 * fjac[i][j+1][k][1][3]
            - tmp1 * njac[i][j+1][k][1][3];
          lhs[i][j][k][CC][1][4] =  tmp2 * fjac[i][j+1][k][1][4]
            - tmp1 * njac[i][j+1][k][1][4];
          lhs[i][j][k][CC][2][0] =  tmp2 * fjac[i][j+1][k][2][0]
            - tmp1 * njac[i][j+1][k][2][0];
          lhs[i][j][k][CC][2][1] =  tmp2 * fjac[i][j+1][k][2][1]
            - tmp1 * njac[i][j+1][k][2][1];
          lhs[i][j][k][CC][2][2] =  tmp2 * fjac[i][j+1][k][2][2]
            - tmp1 * njac[i][j+1][k][2][2]
            - tmp1 * dy3;
          lhs[i][j][k][CC][2][3] =  tmp2 * fjac[i][j+1][k][2][3]
            - tmp1 * njac[i][j+1][k][2][3];
          lhs[i][j][k][CC][2][4] =  tmp2 * fjac[i][j+1][k][2][4]
            - tmp1 * njac[i][j+1][k][2][4];
          lhs[i][j][k][CC][3][0] =  tmp2 * fjac[i][j+1][k][3][0]
            - tmp1 * njac[i][j+1][k][3][0];
          lhs[i][j][k][CC][3][1] =  tmp2 * fjac[i][j+1][k][3][1]
            - tmp1 * njac[i][j+1][k][3][1];
          lhs[i][j][k][CC][3][2] =  tmp2 * fjac[i][j+1][k][3][2]
            - tmp1 * njac[i][j+1][k][3][2];
          lhs[i][j][k][CC][3][3] =  tmp2 * fjac[i][j+1][k][3][3]
            - tmp1 * njac[i][j+1][k][3][3]
            - tmp1 * dy4;
          lhs[i][j][k][CC][3][4] =  tmp2 * fjac[i][j+1][k][3][4]
            - tmp1 * njac[i][j+1][k][3][4];
          lhs[i][j][k][CC][4][0] =  tmp2 * fjac[i][j+1][k][4][0]
            - tmp1 * njac[i][j+1][k][4][0];
          lhs[i][j][k][CC][4][1] =  tmp2 * fjac[i][j+1][k][4][1]
            - tmp1 * njac[i][j+1][k][4][1];
          lhs[i][j][k][CC][4][2] =  tmp2 * fjac[i][j+1][k][4][2]
            - tmp1 * njac[i][j+1][k][4][2];
          lhs[i][j][k][CC][4][3] =  tmp2 * fjac[i][j+1][k][4][3]
            - tmp1 * njac[i][j+1][k][4][3];
          lhs[i][j][k][CC][4][4] =  tmp2 * fjac[i][j+1][k][4][4]
            - tmp1 * njac[i][j+1][k][4][4]
            - tmp1 * dy5;
        }
      }
    }
  } // 结束并行区域
}
// 在文件顶部需要包含 omp.h 头文件
#include <omp.h>
// 假设 pow2 是一个已经定义的宏或函数，例如 #define pow2(x) ((x)*(x))
// 假设所有用到的多维数组和全局变量 (如 u, fjac, njac, lhs, grid_points, c1, c2, dt 等) 
// 都已经在程序的其他地方正确声明和定义。
// 假设 AA, BB, CC 是已定义的枚举或宏。

static void lhsz(void) {
  int i, j, k;

  // 使用一个 parallel 区域来创建一次线程，减少开销
#pragma omp parallel private(i, j, k)
  {
    // ======================================================================
    // 第一个循环：计算 fjac 和 njac 矩阵
    // 使用 omp for 来并行化外层循环。
    // collapse(2) 将 i 和 j 两个循环合并为一个大的迭代空间，有助于更好的负载均衡。
    // 每个线程将分到 i-j 平面的一个子集进行计算。
    // ======================================================================
#pragma omp for collapse(2)
    for (i = 1; i < grid_points[0] - 1; i++) {
      for (j = 1; j < grid_points[1] - 1; j++) {
        for (k = 0; k < grid_points[2]; k++) {
          // 将临时变量在最内层作用域声明，确保它们是每次迭代私有的
          double tmp1, tmp2, tmp3;

          tmp1 = 1.0 / u[i][j][k][0];
          tmp2 = tmp1 * tmp1;
          tmp3 = tmp1 * tmp2;
          
          fjac[i][j][k][0][0] = 0.0;
          fjac[i][j][k][0][1] = 0.0;
          fjac[i][j][k][0][2] = 0.0;
          fjac[i][j][k][0][3] = 1.0;
          fjac[i][j][k][0][4] = 0.0;
          
          fjac[i][j][k][1][0] = - (u[i][j][k][1] * u[i][j][k][3]) * tmp2;
          fjac[i][j][k][1][1] = u[i][j][k][3] * tmp1;
          fjac[i][j][k][1][2] = 0.0;
          fjac[i][j][k][1][3] = u[i][j][k][1] * tmp1;
          fjac[i][j][k][1][4] = 0.0;
          
          fjac[i][j][k][2][0] = - (u[i][j][k][2] * u[i][j][k][3]) * tmp2;
          fjac[i][j][k][2][1] = 0.0;
          fjac[i][j][k][2][2] = u[i][j][k][3] * tmp1;
          fjac[i][j][k][2][3] = u[i][j][k][2] * tmp1;
          fjac[i][j][k][2][4] = 0.0;
          
          fjac[i][j][k][3][0] = - (u[i][j][k][3] * u[i][j][k][3] * tmp2) + 0.50 * c2 * ((u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2);
          fjac[i][j][k][3][1] = - c2 * u[i][j][k][1] * tmp1;
          fjac[i][j][k][3][2] = - c2 * u[i][j][k][2] * tmp1;
          fjac[i][j][k][3][3] = (2.0 - c2) * u[i][j][k][3] * tmp1;
          fjac[i][j][k][3][4] = c2;
          
          fjac[i][j][k][4][0] = (c2 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2 - c1 * (u[i][j][k][4] * tmp1)) * (u[i][j][k][3] * tmp1);
          fjac[i][j][k][4][1] = - c2 * (u[i][j][k][1] * u[i][j][k][3]) * tmp2;
          fjac[i][j][k][4][2] = - c2 * (u[i][j][k][2] * u[i][j][k][3]) * tmp2;
          fjac[i][j][k][4][3] = c1 * (u[i][j][k][4] * tmp1) - 0.50 * c2 * ((u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + 3.0 * u[i][j][k][3] * u[i][j][k][3]) * tmp2);
          fjac[i][j][k][4][4] = c1 * u[i][j][k][3] * tmp1;
          
          njac[i][j][k][0][0] = 0.0;
          njac[i][j][k][0][1] = 0.0;
          njac[i][j][k][0][2] = 0.0;
          njac[i][j][k][0][3] = 0.0;
          njac[i][j][k][0][4] = 0.0;
          
          njac[i][j][k][1][0] = - c3c4 * tmp2 * u[i][j][k][1];
          njac[i][j][k][1][1] = c3c4 * tmp1;
          njac[i][j][k][1][2] = 0.0;
          njac[i][j][k][1][3] = 0.0;
          njac[i][j][k][1][4] = 0.0;
          
          njac[i][j][k][2][0] = - c3c4 * tmp2 * u[i][j][k][2];
          njac[i][j][k][2][1] = 0.0;
          njac[i][j][k][2][2] = c3c4 * tmp1;
          njac[i][j][k][2][3] = 0.0;
          njac[i][j][k][2][4] = 0.0;
          
          njac[i][j][k][3][0] = - con43 * c3c4 * tmp2 * u[i][j][k][3];
          njac[i][j][k][3][1] = 0.0;
          njac[i][j][k][3][2] = 0.0;
          njac[i][j][k][3][3] = con43 * c3 * c4 * tmp1;
          njac[i][j][k][3][4] = 0.0;
          
          njac[i][j][k][4][0] = - (c3c4 - c1345) * tmp3 * (pow2(u[i][j][k][1])) - (c3c4 - c1345) * tmp3 * (pow2(u[i][j][k][2])) - (con43 * c3c4 - c1345) * tmp3 * (pow2(u[i][j][k][3])) - c1345 * tmp2 * u[i][j][k][4];
          njac[i][j][k][4][1] = (c3c4 - c1345) * tmp2 * u[i][j][k][1];
          njac[i][j][k][4][2] = (c3c4 - c1345) * tmp2 * u[i][j][k][2];
          njac[i][j][k][4][3] = (con43 * c3c4 - c1345) * tmp2 * u[i][j][k][3];
          njac[i][j][k][4][4] = (c1345) * tmp1;
        }
      }
    }

    // ======================================================================
    // 此处有一个隐式屏障 (implicit barrier)。
    // 所有线程在完成上面的循环后会在此等待，直到所有线程都完成。
    // 这确保了在进入下一个循环之前，fjac 和 njac 的所有值都已计算完毕。
    // ======================================================================

    // ======================================================================
    // 第二个循环：计算 lhs 矩阵
    // 同样，我们使用 omp for 和 collapse(2) 来并行化这个循环。
    // ======================================================================
#pragma omp for collapse(2)
    for (i = 1; i < grid_points[0] - 1; i++) {
      for (j = 1; j < grid_points[1] - 1; j++) {
        for (k = 1; k < grid_points[2] - 1; k++) {
          // 将临时变量在最内层作用域声明，确保它们是每次迭代私有的
          double tmp1, tmp2;

          tmp1 = dt * tz1;
          tmp2 = dt * tz2;
          
          lhs[i][j][k][AA][0][0] = - tmp2 * fjac[i][j][k-1][0][0] - tmp1 * njac[i][j][k-1][0][0] - tmp1 * dz1;
          lhs[i][j][k][AA][0][1] = - tmp2 * fjac[i][j][k-1][0][1] - tmp1 * njac[i][j][k-1][0][1];
          lhs[i][j][k][AA][0][2] = - tmp2 * fjac[i][j][k-1][0][2] - tmp1 * njac[i][j][k-1][0][2];
          lhs[i][j][k][AA][0][3] = - tmp2 * fjac[i][j][k-1][0][3] - tmp1 * njac[i][j][k-1][0][3];
          lhs[i][j][k][AA][0][4] = - tmp2 * fjac[i][j][k-1][0][4] - tmp1 * njac[i][j][k-1][0][4];
          
          lhs[i][j][k][AA][1][0] = - tmp2 * fjac[i][j][k-1][1][0] - tmp1 * njac[i][j][k-1][1][0];
          lhs[i][j][k][AA][1][1] = - tmp2 * fjac[i][j][k-1][1][1] - tmp1 * njac[i][j][k-1][1][1] - tmp1 * dz2;
          lhs[i][j][k][AA][1][2] = - tmp2 * fjac[i][j][k-1][1][2] - tmp1 * njac[i][j][k-1][1][2];
          lhs[i][j][k][AA][1][3] = - tmp2 * fjac[i][j][k-1][1][3] - tmp1 * njac[i][j][k-1][1][3];
          lhs[i][j][k][AA][1][4] = - tmp2 * fjac[i][j][k-1][1][4] - tmp1 * njac[i][j][k-1][1][4];
          
          lhs[i][j][k][AA][2][0] = - tmp2 * fjac[i][j][k-1][2][0] - tmp1 * njac[i][j][k-1][2][0];
          lhs[i][j][k][AA][2][1] = - tmp2 * fjac[i][j][k-1][2][1] - tmp1 * njac[i][j][k-1][2][1];
          lhs[i][j][k][AA][2][2] = - tmp2 * fjac[i][j][k-1][2][2] - tmp1 * njac[i][j][k-1][2][2] - tmp1 * dz3;
          lhs[i][j][k][AA][2][3] = - tmp2 * fjac[i][j][k-1][2][3] - tmp1 * njac[i][j][k-1][2][3];
          lhs[i][j][k][AA][2][4] = - tmp2 * fjac[i][j][k-1][2][4] - tmp1 * njac[i][j][k-1][2][4];
          
          lhs[i][j][k][AA][3][0] = - tmp2 * fjac[i][j][k-1][3][0] - tmp1 * njac[i][j][k-1][3][0];
          lhs[i][j][k][AA][3][1] = - tmp2 * fjac[i][j][k-1][3][1] - tmp1 * njac[i][j][k-1][3][1];
          lhs[i][j][k][AA][3][2] = - tmp2 * fjac[i][j][k-1][3][2] - tmp1 * njac[i][j][k-1][3][2];
          lhs[i][j][k][AA][3][3] = - tmp2 * fjac[i][j][k-1][3][3] - tmp1 * njac[i][j][k-1][3][3] - tmp1 * dz4;
          lhs[i][j][k][AA][3][4] = - tmp2 * fjac[i][j][k-1][3][4] - tmp1 * njac[i][j][k-1][3][4];
          
          lhs[i][j][k][AA][4][0] = - tmp2 * fjac[i][j][k-1][4][0] - tmp1 * njac[i][j][k-1][4][0];
          lhs[i][j][k][AA][4][1] = - tmp2 * fjac[i][j][k-1][4][1] - tmp1 * njac[i][j][k-1][4][1];
          lhs[i][j][k][AA][4][2] = - tmp2 * fjac[i][j][k-1][4][2] - tmp1 * njac[i][j][k-1][4][2];
          lhs[i][j][k][AA][4][3] = - tmp2 * fjac[i][j][k-1][4][3] - tmp1 * njac[i][j][k-1][4][3];
          lhs[i][j][k][AA][4][4] = - tmp2 * fjac[i][j][k-1][4][4] - tmp1 * njac[i][j][k-1][4][4] - tmp1 * dz5;
          
          lhs[i][j][k][BB][0][0] = 1.0 + tmp1 * 2.0 * njac[i][j][k][0][0] + tmp1 * 2.0 * dz1;
          lhs[i][j][k][BB][0][1] = tmp1 * 2.0 * njac[i][j][k][0][1];
          lhs[i][j][k][BB][0][2] = tmp1 * 2.0 * njac[i][j][k][0][2];
          lhs[i][j][k][BB][0][3] = tmp1 * 2.0 * njac[i][j][k][0][3];
          lhs[i][j][k][BB][0][4] = tmp1 * 2.0 * njac[i][j][k][0][4];
          
          lhs[i][j][k][BB][1][0] = tmp1 * 2.0 * njac[i][j][k][1][0];
          lhs[i][j][k][BB][1][1] = 1.0 + tmp1 * 2.0 * njac[i][j][k][1][1] + tmp1 * 2.0 * dz2;
          lhs[i][j][k][BB][1][2] = tmp1 * 2.0 * njac[i][j][k][1][2];
          lhs[i][j][k][BB][1][3] = tmp1 * 2.0 * njac[i][j][k][1][3];
          lhs[i][j][k][BB][1][4] = tmp1 * 2.0 * njac[i][j][k][1][4];
          
          lhs[i][j][k][BB][2][0] = tmp1 * 2.0 * njac[i][j][k][2][0];
          lhs[i][j][k][BB][2][1] = tmp1 * 2.0 * njac[i][j][k][2][1];
          lhs[i][j][k][BB][2][2] = 1.0 + tmp1 * 2.0 * njac[i][j][k][2][2] + tmp1 * 2.0 * dz3;
          lhs[i][j][k][BB][2][3] = tmp1 * 2.0 * njac[i][j][k][2][3];
          lhs[i][j][k][BB][2][4] = tmp1 * 2.0 * njac[i][j][k][2][4];
          
          lhs[i][j][k][BB][3][0] = tmp1 * 2.0 * njac[i][j][k][3][0];
          lhs[i][j][k][BB][3][1] = tmp1 * 2.0 * njac[i][j][k][3][1];
          lhs[i][j][k][BB][3][2] = tmp1 * 2.0 * njac[i][j][k][3][2];
          lhs[i][j][k][BB][3][3] = 1.0 + tmp1 * 2.0 * njac[i][j][k][3][3] + tmp1 * 2.0 * dz4;
          lhs[i][j][k][BB][3][4] = tmp1 * 2.0 * njac[i][j][k][3][4];
          
          lhs[i][j][k][BB][4][0] = tmp1 * 2.0 * njac[i][j][k][4][0];
          lhs[i][j][k][BB][4][1] = tmp1 * 2.0 * njac[i][j][k][4][1];
          lhs[i][j][k][BB][4][2] = tmp1 * 2.0 * njac[i][j][k][4][2];
          lhs[i][j][k][BB][4][3] = tmp1 * 2.0 * njac[i][j][k][4][3];
          lhs[i][j][k][BB][4][4] = 1.0 + tmp1 * 2.0 * njac[i][j][k][4][4] + tmp1 * 2.0 * dz5;
          
          lhs[i][j][k][CC][0][0] = tmp2 * fjac[i][j][k+1][0][0] - tmp1 * njac[i][j][k+1][0][0] - tmp1 * dz1;
          lhs[i][j][k][CC][0][1] = tmp2 * fjac[i][j][k+1][0][1] - tmp1 * njac[i][j][k+1][0][1];
          lhs[i][j][k][CC][0][2] = tmp2 * fjac[i][j][k+1][0][2] - tmp1 * njac[i][j][k+1][0][2];
          lhs[i][j][k][CC][0][3] = tmp2 * fjac[i][j][k+1][0][3] - tmp1 * njac[i][j][k+1][0][3];
          lhs[i][j][k][CC][0][4] = tmp2 * fjac[i][j][k+1][0][4] - tmp1 * njac[i][j][k+1][0][4];
          
          lhs[i][j][k][CC][1][0] = tmp2 * fjac[i][j][k+1][1][0] - tmp1 * njac[i][j][k+1][1][0];
          lhs[i][j][k][CC][1][1] = tmp2 * fjac[i][j][k+1][1][1] - tmp1 * njac[i][j][k+1][1][1] - tmp1 * dz2;
          lhs[i][j][k][CC][1][2] = tmp2 * fjac[i][j][k+1][1][2] - tmp1 * njac[i][j][k+1][1][2];
          lhs[i][j][k][CC][1][3] = tmp2 * fjac[i][j][k+1][1][3] - tmp1 * njac[i][j][k+1][1][3];
          lhs[i][j][k][CC][1][4] = tmp2 * fjac[i][j][k+1][1][4] - tmp1 * njac[i][j][k+1][1][4];
          
          lhs[i][j][k][CC][2][0] = tmp2 * fjac[i][j][k+1][2][0] - tmp1 * njac[i][j][k+1][2][0];
          lhs[i][j][k][CC][2][1] = tmp2 * fjac[i][j][k+1][2][1] - tmp1 * njac[i][j][k+1][2][1];
          lhs[i][j][k][CC][2][2] = tmp2 * fjac[i][j][k+1][2][2] - tmp1 * njac[i][j][k+1][2][2] - tmp1 * dz3;
          lhs[i][j][k][CC][2][3] = tmp2 * fjac[i][j][k+1][2][3] - tmp1 * njac[i][j][k+1][2][3];
          lhs[i][j][k][CC][2][4] = tmp2 * fjac[i][j][k+1][2][4] - tmp1 * njac[i][j][k+1][2][4];
          
          lhs[i][j][k][CC][3][0] = tmp2 * fjac[i][j][k+1][3][0] - tmp1 * njac[i][j][k+1][3][0];
          lhs[i][j][k][CC][3][1] = tmp2 * fjac[i][j][k+1][3][1] - tmp1 * njac[i][j][k+1][3][1];
          lhs[i][j][k][CC][3][2] = tmp2 * fjac[i][j][k+1][3][2] - tmp1 * njac[i][j][k+1][3][2];
          lhs[i][j][k][CC][3][3] = tmp2 * fjac[i][j][k+1][3][3] - tmp1 * njac[i][j][k+1][3][3] - tmp1 * dz4;
          lhs[i][j][k][CC][3][4] = tmp2 * fjac[i][j][k+1][3][4] - tmp1 * njac[i][j][k+1][3][4];
          
          lhs[i][j][k][CC][4][0] = tmp2 * fjac[i][j][k+1][4][0] - tmp1 * njac[i][j][k+1][4][0];
          lhs[i][j][k][CC][4][1] = tmp2 * fjac[i][j][k+1][4][1] - tmp1 * njac[i][j][k+1][4][1];
          lhs[i][j][k][CC][4][2] = tmp2 * fjac[i][j][k+1][4][2] - tmp1 * njac[i][j][k+1][4][2];
          lhs[i][j][k][CC][4][3] = tmp2 * fjac[i][j][k+1][4][3] - tmp1 * njac[i][j][k+1][4][3];
          lhs[i][j][k][CC][4][4] = tmp2 * fjac[i][j][k+1][4][4] - tmp1 * njac[i][j][k+1][4][4] - tmp1 * dz5;
        }
      }
    }

  } // end of parallel region
}
#include <omp.h> // 使用 OpenMP 魔法，需要先告诉电脑

static void compute_rhs(void) {
  int i, j, k, m;
  // 这些是临时用的小工具，每个小厨师都需要自己的一套
  double rho_inv, uijk, up1, um1, vijk, vp1, vm1, wijk, wp1, wm1;

  // 魔法开始！把整个函数变成一个大的团队合作区
  // private(...)：括号里的是每个小厨师自己独有的工具
  #pragma omp parallel private(i, j, k, m, rho_inv, uijk, up1, um1, vijk, vp1, vm1, wijk, wp1, wm1)
  {
    // --- 第1步：准备工作 ---
    // collapse(3)：把三层循环合并成一个大任务分给大家
    #pragma omp for collapse(3)
    for (i = 0; i < grid_points[0]; i++) {
      for (j = 0; j < grid_points[1]; j++) {
        for (k = 0; k < grid_points[2]; k++) {
          rho_inv = 1.0/u[i][j][k][0];
          rho_i[i][j][k] = rho_inv;
          us[i][j][k] = u[i][j][k][1] * rho_inv;
          vs[i][j][k] = u[i][j][k][2] * rho_inv;
          ws[i][j][k] = u[i][j][k][3] * rho_inv;
          square[i][j][k] = 0.5 * (u[i][j][k][1]*u[i][j][k][1] + 
                   u[i][j][k][2]*u[i][j][k][2] +
                   u[i][j][k][3]*u[i][j][k][3] ) * rho_inv;
          qs[i][j][k] = square[i][j][k] * rho_inv;
        }
      }
    }

    // --- 第2步：把一个盘子里的菜（forcing）倒到另一个盘子里（rhs） ---
    #pragma omp for collapse(3)
    for (i = 0; i < grid_points[0]; i++) {
      for (j = 0; j < grid_points[1]; j++) {
        for (k = 0; k < grid_points[2]; k++) {
          for (m = 0; m < 5; m++) {
            rhs[i][j][k][m] = forcing[i][j][k][m];
          }
        }
      }
    }

    // --- 第3步：X方向（可以想成是左右方向）的计算 ---
    #pragma omp for collapse(3)
    for (i = 1; i < grid_points[0]-1; i++) {
      for (j = 1; j < grid_points[1]-1; j++) {
        for (k = 1; k < grid_points[2]-1; k++) {
          uijk = us[i][j][k];
          up1  = us[i+1][j][k];
          um1  = us[i-1][j][k];
          rhs[i][j][k][0] = rhs[i][j][k][0] + dx1tx1 * 
            (u[i+1][j][k][0] - 2.0*u[i][j][k][0] + 
             u[i-1][j][k][0]) -
            tx2 * (u[i+1][j][k][1] - u[i-1][j][k][1]);
          // ... (为了简洁，省略了中间一长串一模一样的计算)
          rhs[i][j][k][1] = rhs[i][j][k][1] + dx2tx1 * (u[i+1][j][k][1] - 2.0*u[i][j][k][1] + u[i-1][j][k][1]) + xxcon2*con43 * (up1 - 2.0*uijk + um1) - tx2 * (u[i+1][j][k][1]*up1 - u[i-1][j][k][1]*um1 + (u[i+1][j][k][4]- square[i+1][j][k]- u[i-1][j][k][4]+ square[i-1][j][k])* c2);
          rhs[i][j][k][2] = rhs[i][j][k][2] + dx3tx1 * (u[i+1][j][k][2] - 2.0*u[i][j][k][2] + u[i-1][j][k][2]) + xxcon2 * (vs[i+1][j][k] - 2.0*vs[i][j][k] + vs[i-1][j][k]) - tx2 * (u[i+1][j][k][2]*up1 - u[i-1][j][k][2]*um1);
          rhs[i][j][k][3] = rhs[i][j][k][3] + dx4tx1 * (u[i+1][j][k][3] - 2.0*u[i][j][k][3] + u[i-1][j][k][3]) + xxcon2 * (ws[i+1][j][k] - 2.0*ws[i][j][k] + ws[i-1][j][k]) - tx2 * (u[i+1][j][k][3]*up1 - u[i-1][j][k][3]*um1);
          rhs[i][j][k][4] = rhs[i][j][k][4] + dx5tx1 * (u[i+1][j][k][4] - 2.0*u[i][j][k][4] + u[i-1][j][k][4]) + xxcon3 * (qs[i+1][j][k] - 2.0*qs[i][j][k] + qs[i-1][j][k]) + xxcon4 * (up1*up1 - 2.0*uijk*uijk + um1*um1) + xxcon5 * (u[i+1][j][k][4]*rho_i[i+1][j][k] - 2.0*u[i][j][k][4]*rho_i[i][j][k] + u[i-1][j][k][4]*rho_i[i-1][j][k]) - tx2 * ( (c1*u[i+1][j][k][4] - c2*square[i+1][j][k])*up1 - (c1*u[i-1][j][k][4] - c2*square[i-1][j][k])*um1 );
        }
      }
    }
    
    // --- X方向边界处理 ---
    // 这里的 i 是固定的，所以我们在 j 和 k 上分配任务
    #pragma omp for collapse(2)
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
          i = 1;
          for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m]- dssp * ( 5.0*u[i][j][k][m] - 4.0*u[i+1][j][k][m] + u[i+2][j][k][m]); }
      }
    }
    #pragma omp for collapse(2)
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
          i = 2;
          for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (-4.0*u[i-1][j][k][m] + 6.0*u[i][j][k][m] - 4.0*u[i+1][j][k][m] + u[i+2][j][k][m]); }
      }
    }
    #pragma omp for collapse(2)
    for (i = 3; i < grid_points[0]-3; i++) {
      for (j = 1; j < grid_points[1]-1; j++) {
          for (k = 1; k < grid_points[2]-1; k++) {
              for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * ( u[i-2][j][k][m] - 4.0*u[i-1][j][k][m] + 6.0*u[i][j][k][m] - 4.0*u[i+1][j][k][m] + u[i+2][j][k][m] ); }
          }
      }
    }
    #pragma omp for collapse(2)
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
          i = grid_points[0]-3;
          for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * ( u[i-2][j][k][m] - 4.0*u[i-1][j][k][m] + 6.0*u[i][j][k][m] - 4.0*u[i+1][j][k][m] ); }
      }
    }
    #pragma omp for collapse(2)
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
          i = grid_points[0]-2;
          for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * ( u[i-2][j][k][m] - 4.*u[i-1][j][k][m] + 5.0*u[i][j][k][m] ); }
      }
    }

    // --- 第4步：Y方向（可以想成是前后方向）的计算 ---
    // 这里的逻辑和X方向几乎一样，只是操作的方向变了
    #pragma omp for collapse(3)
    for (i = 1; i < grid_points[0]-1; i++) {
      for (j = 1; j < grid_points[1]-1; j++) {
        for (k = 1; k < grid_points[2]-1; k++) {
          vijk = vs[i][j][k];
          vp1  = vs[i][j+1][k];
          vm1  = vs[i][j-1][k];
          // ... (省略了中间一长串一模一样的计算)
          rhs[i][j][k][0] = rhs[i][j][k][0] + dy1ty1 * (u[i][j+1][k][0] - 2.0*u[i][j][k][0] + u[i][j-1][k][0]) - ty2 * (u[i][j+1][k][2] - u[i][j-1][k][2]);
          rhs[i][j][k][1] = rhs[i][j][k][1] + dy2ty1 * (u[i][j+1][k][1] - 2.0*u[i][j][k][1] + u[i][j-1][k][1]) + yycon2 * (us[i][j+1][k] - 2.0*us[i][j][k] + us[i][j-1][k]) - ty2 * (u[i][j+1][k][1]*vp1 - u[i][j-1][k][1]*vm1);
          rhs[i][j][k][2] = rhs[i][j][k][2] + dy3ty1 * (u[i][j+1][k][2] - 2.0*u[i][j][k][2] + u[i][j-1][k][2]) + yycon2*con43 * (vp1 - 2.0*vijk + vm1) - ty2 * (u[i][j+1][k][2]*vp1 - u[i][j-1][k][2]*vm1 + (u[i][j+1][k][4] - square[i][j+1][k] - u[i][j-1][k][4] + square[i][j-1][k]) *c2);
          rhs[i][j][k][3] = rhs[i][j][k][3] + dy4ty1 * (u[i][j+1][k][3] - 2.0*u[i][j][k][3] + u[i][j-1][k][3]) + yycon2 * (ws[i][j+1][k] - 2.0*ws[i][j][k] + ws[i][j-1][k]) - ty2 * (u[i][j+1][k][3]*vp1 - u[i][j-1][k][3]*vm1);
          rhs[i][j][k][4] = rhs[i][j][k][4] + dy5ty1 * (u[i][j+1][k][4] - 2.0*u[i][j][k][4] + u[i][j-1][k][4]) + yycon3 * (qs[i][j+1][k] - 2.0*qs[i][j][k] + qs[i][j-1][k]) + yycon4 * (vp1*vp1 - 2.0*vijk*vijk + vm1*vm1) + yycon5 * (u[i][j+1][k][4]*rho_i[i][j+1][k] - 2.0*u[i][j][k][4]*rho_i[i][j][k] + u[i][j-1][k][4]*rho_i[i][j-1][k]) - ty2 * ((c1*u[i][j+1][k][4] - c2*square[i][j+1][k]) * vp1 - (c1*u[i][j-1][k][4] - c2*square[i][j-1][k]) * vm1);
        }
      }
    }

    // --- Y方向边界处理 ---
    #pragma omp for collapse(2)
    for (i = 1; i < grid_points[0]-1; i++) {
        for (k = 1; k < grid_points[2]-1; k++) {
            j = 1;
            for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m]- dssp * ( 5.0*u[i][j][k][m] - 4.0*u[i][j+1][k][m] + u[i][j+2][k][m]); }
        }
    }
    #pragma omp for collapse(2)
    for (i = 1; i < grid_points[0]-1; i++) {
        for (k = 1; k < grid_points[2]-1; k++) {
            j = 2;
            for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (-4.0*u[i][j-1][k][m] + 6.0*u[i][j][k][m] - 4.0*u[i][j+1][k][m] + u[i][j+2][k][m]); }
        }
    }
    #pragma omp for collapse(3)
    for (i = 1; i < grid_points[0]-1; i++) {
        for (j = 3; j < grid_points[1]-3; j++) {
            for (k = 1; k < grid_points[2]-1; k++) {
                for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * ( u[i][j-2][k][m] - 4.0*u[i][j-1][k][m] + 6.0*u[i][j][k][m] - 4.0*u[i][j+1][k][m] + u[i][j+2][k][m] ); }
            }
        }
    }
    #pragma omp for collapse(2)
    for (i = 1; i < grid_points[0]-1; i++) {
        for (k = 1; k < grid_points[2]-1; k++) {
            j = grid_points[1]-3;
            for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * ( u[i][j-2][k][m] - 4.0*u[i][j-1][k][m] + 6.0*u[i][j][k][m] - 4.0*u[i][j+1][k][m] ); }
        }
    }
    #pragma omp for collapse(2)
    for (i = 1; i < grid_points[0]-1; i++) {
        for (k = 1; k < grid_points[2]-1; k++) {
            j = grid_points[1]-2;
            for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * ( u[i][j-2][k][m] - 4.*u[i][j-1][k][m] + 5.*u[i][j][k][m] ); }
        }
    }

    // --- 第5步：Z方向（可以想成是上下方向）的计算 ---
    #pragma omp for collapse(3)
    for (i = 1; i < grid_points[0]-1; i++) {
      for (j = 1; j < grid_points[1]-1; j++) {
        for (k = 1; k < grid_points[2]-1; k++) {
          wijk = ws[i][j][k];
          wp1  = ws[i][j][k+1];
          wm1  = ws[i][j][k-1];
          // ... (省略了中间一长串一模一样的计算)
          rhs[i][j][k][0] = rhs[i][j][k][0] + dz1tz1 * (u[i][j][k+1][0] - 2.0*u[i][j][k][0] + u[i][j][k-1][0]) - tz2 * (u[i][j][k+1][3] - u[i][j][k-1][3]);
          rhs[i][j][k][1] = rhs[i][j][k][1] + dz2tz1 * (u[i][j][k+1][1] - 2.0*u[i][j][k][1] + u[i][j][k-1][1]) + zzcon2 * (us[i][j][k+1] - 2.0*us[i][j][k] + us[i][j][k-1]) - tz2 * (u[i][j][k+1][1]*wp1 - u[i][j][k-1][1]*wm1);
          rhs[i][j][k][2] = rhs[i][j][k][2] + dz3tz1 * (u[i][j][k+1][2] - 2.0*u[i][j][k][2] + u[i][j][k-1][2]) + zzcon2 * (vs[i][j][k+1] - 2.0*vs[i][j][k] + vs[i][j][k-1]) - tz2 * (u[i][j][k+1][2]*wp1 - u[i][j][k-1][2]*wm1);
          rhs[i][j][k][3] = rhs[i][j][k][3] + dz4tz1 * (u[i][j][k+1][3] - 2.0*u[i][j][k][3] + u[i][j][k-1][3]) + zzcon2*con43 * (wp1 - 2.0*wijk + wm1) - tz2 * (u[i][j][k+1][3]*wp1 - u[i][j][k-1][3]*wm1 + (u[i][j][k+1][4] - square[i][j][k+1] - u[i][j][k-1][4] + square[i][j][k-1]) *c2);
          rhs[i][j][k][4] = rhs[i][j][k][4] + dz5tz1 * (u[i][j][k+1][4] - 2.0*u[i][j][k][4] + u[i][j][k-1][4]) + zzcon3 * (qs[i][j][k+1] - 2.0*qs[i][j][k] + qs[i][j][k-1]) + zzcon4 * (wp1*wp1 - 2.0*wijk*wijk + wm1*wm1) + zzcon5 * (u[i][j][k+1][4]*rho_i[i][j][k+1] - 2.0*u[i][j][k][4]*rho_i[i][j][k] + u[i][j][k-1][4]*rho_i[i][j][k-1]) - tz2 * ( (c1*u[i][j][k+1][4] - c2*square[i][j][k+1])*wp1 - (c1*u[i][j][k-1][4] - c2*square[i][j][k-1])*wm1);
        }
      }
    }
    
    // --- Z方向边界处理 ---
    #pragma omp for collapse(2)
    for (i = 1; i < grid_points[0]-1; i++) {
        for (j = 1; j < grid_points[1]-1; j++) {
            k = 1;
            for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m]- dssp * ( 5.0*u[i][j][k][m] - 4.0*u[i][j][k+1][m] + u[i][j][k+2][m]); }
        }
    }
    #pragma omp for collapse(2)
    for (i = 1; i < grid_points[0]-1; i++) {
        for (j = 1; j < grid_points[1]-1; j++) {
            k = 2;
            for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (-4.0*u[i][j][k-1][m] + 6.0*u[i][j][k][m] - 4.0*u[i][j][k+1][m] + u[i][j][k+2][m]); }
        }
    }
    #pragma omp for collapse(3)
    for (i = 1; i < grid_points[0]-1; i++) {
        for (j = 1; j < grid_points[1]-1; j++) {
            for (k = 3; k < grid_points[2]-3; k++) {
                for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * ( u[i][j][k-2][m] - 4.0*u[i][j][k-1][m] + 6.0*u[i][j][k][m] - 4.0*u[i][j][k+1][m] + u[i][j][k+2][m] ); }
            }
        }
    }
    #pragma omp for collapse(2)
    for (i = 1; i < grid_points[0]-1; i++) {
        for (j = 1; j < grid_points[1]-1; j++) {
            k = grid_points[2]-3;
            for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * ( u[i][j][k-2][m] - 4.0*u[i][j][k-1][m] + 6.0*u[i][j][k][m] - 4.0*u[i][j][k+1][m] ); }
        }
    }
    #pragma omp for collapse(2)
    for (i = 1; i < grid_points[0]-1; i++) {
        for (j = 1; j < grid_points[1]-1; j++) {
            k = grid_points[2]-2;
            for (m = 0; m < 5; m++) { rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * ( u[i][j][k-2][m] - 4.0*u[i][j][k-1][m] + 5.0*u[i][j][k][m] ); }
        }
    }
    
    // --- 第6步：最后的加工，给所有菜都乘上一个系数 ---
    // 这里我们把 j 和 k 合并起来分工
    #pragma omp for collapse(2)
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
        for (m = 0; m < 5; m++) {
          for (i = 1; i < grid_points[0]-1; i++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] * dt;
          }
        }
      }
    }

  } // 团队合作区结束，小厨师们可以休息了
}
static void set_constants(void) {
  ce[0][0]  = 2.0;
  ce[0][1]  = 0.0;
  ce[0][2]  = 0.0;
  ce[0][3]  = 4.0;
  ce[0][4]  = 5.0;
  ce[0][5]  = 3.0;
  ce[0][6]  = 0.5;
  ce[0][7]  = 0.02;
  ce[0][8]  = 0.01;
  ce[0][9]  = 0.03;
  ce[0][10] = 0.5;
  ce[0][11] = 0.4;
  ce[0][12] = 0.3;
  ce[1][0]  = 1.0;
  ce[1][1]  = 0.0;
  ce[1][2]  = 0.0;
  ce[1][3]  = 0.0;
  ce[1][4]  = 1.0;
  ce[1][5]  = 2.0;
  ce[1][6]  = 3.0;
  ce[1][7]  = 0.01;
  ce[1][8]  = 0.03;
  ce[1][9]  = 0.02;
  ce[1][10] = 0.4;
  ce[1][11] = 0.3;
  ce[1][12] = 0.5;
  ce[2][0]  = 2.0;
  ce[2][1]  = 2.0;
  ce[2][2]  = 0.0;
  ce[2][3]  = 0.0;
  ce[2][4]  = 0.0;
  ce[2][5]  = 2.0;
  ce[2][6]  = 3.0;
  ce[2][7]  = 0.04;
  ce[2][8]  = 0.03;
  ce[2][9]  = 0.05;
  ce[2][10] = 0.3;
  ce[2][11] = 0.5;
  ce[2][12] = 0.4;
  ce[3][0]  = 2.0;
  ce[3][1]  = 2.0;
  ce[3][2]  = 0.0;
  ce[3][3]  = 0.0;
  ce[3][4]  = 0.0;
  ce[3][5]  = 2.0;
  ce[3][6]  = 3.0;
  ce[3][7]  = 0.03;
  ce[3][8]  = 0.05;
  ce[3][9]  = 0.04;
  ce[3][10] = 0.2;
  ce[3][11] = 0.1;
  ce[3][12] = 0.3;
  ce[4][0]  = 5.0;
  ce[4][1]  = 4.0;
  ce[4][2]  = 3.0;
  ce[4][3]  = 2.0;
  ce[4][4]  = 0.1;
  ce[4][5]  = 0.4;
  ce[4][6]  = 0.3;
  ce[4][7]  = 0.05;
  ce[4][8]  = 0.04;
  ce[4][9]  = 0.03;
  ce[4][10] = 0.1;
  ce[4][11] = 0.3;
  ce[4][12] = 0.2;
  c1 = 1.4;
  c2 = 0.4;
  c3 = 0.1;
  c4 = 1.0;
  c5 = 1.4;
  dnxm1 = 1.0 / (double)(grid_points[0]-1);
  dnym1 = 1.0 / (double)(grid_points[1]-1);
  dnzm1 = 1.0 / (double)(grid_points[2]-1);
  c1c2 = c1 * c2;
  c1c5 = c1 * c5;
  c3c4 = c3 * c4;
  c1345 = c1c5 * c3c4;
  conz1 = (1.0-c1c5);
  tx1 = 1.0 / (dnxm1 * dnxm1);
  tx2 = 1.0 / (2.0 * dnxm1);
  tx3 = 1.0 / dnxm1;
  ty1 = 1.0 / (dnym1 * dnym1);
  ty2 = 1.0 / (2.0 * dnym1);
  ty3 = 1.0 / dnym1;
  tz1 = 1.0 / (dnzm1 * dnzm1);
  tz2 = 1.0 / (2.0 * dnzm1);
  tz3 = 1.0 / dnzm1;
  dx1 = 0.75;
  dx2 = 0.75;
  dx3 = 0.75;
  dx4 = 0.75;
  dx5 = 0.75;
  dy1 = 0.75;
  dy2 = 0.75;
  dy3 = 0.75;
  dy4 = 0.75;
  dy5 = 0.75;
  dz1 = 1.0;
  dz2 = 1.0;
  dz3 = 1.0;
  dz4 = 1.0;
  dz5 = 1.0;
  dxmax = max(dx3, dx4);
  dymax = max(dy2, dy4);
  dzmax = max(dz2, dz3);
  dssp = 0.25 * max(dx1, max(dy1, dz1) );
  c4dssp = 4.0 * dssp;
  c5dssp = 5.0 * dssp;
  dttx1 = dt*tx1;
  dttx2 = dt*tx2;
  dtty1 = dt*ty1;
  dtty2 = dt*ty2;
  dttz1 = dt*tz1;
  dttz2 = dt*tz2;
  c2dttx1 = 2.0*dttx1;
  c2dtty1 = 2.0*dtty1;
  c2dttz1 = 2.0*dttz1;
  dtdssp = dt*dssp;
  comz1  = dtdssp;
  comz4  = 4.0*dtdssp;
  comz5  = 5.0*dtdssp;
  comz6  = 6.0*dtdssp;
  c3c4tx3 = c3c4*tx3;
  c3c4ty3 = c3c4*ty3;
  c3c4tz3 = c3c4*tz3;
  dx1tx1 = dx1*tx1;
  dx2tx1 = dx2*tx1;
  dx3tx1 = dx3*tx1;
  dx4tx1 = dx4*tx1;
  dx5tx1 = dx5*tx1;
  dy1ty1 = dy1*ty1;
  dy2ty1 = dy2*ty1;
  dy3ty1 = dy3*ty1;
  dy4ty1 = dy4*ty1;
  dy5ty1 = dy5*ty1;
  dz1tz1 = dz1*tz1;
  dz2tz1 = dz2*tz1;
  dz3tz1 = dz3*tz1;
  dz4tz1 = dz4*tz1;
  dz5tz1 = dz5*tz1;
  c2iv  = 2.5;
  con43 = 4.0/3.0;
  con16 = 1.0/6.0;
  xxcon1 = c3c4tx3*con43*tx3;
  xxcon2 = c3c4tx3*tx3;
  xxcon3 = c3c4tx3*conz1*tx3;
  xxcon4 = c3c4tx3*con16*tx3;
  xxcon5 = c3c4tx3*c1c5*tx3;
  yycon1 = c3c4ty3*con43*ty3;
  yycon2 = c3c4ty3*ty3;
  yycon3 = c3c4ty3*conz1*ty3;
  yycon4 = c3c4ty3*con16*ty3;
  yycon5 = c3c4ty3*c1c5*ty3;
  zzcon1 = c3c4tz3*con43*tz3;
  zzcon2 = c3c4tz3*tz3;
  zzcon3 = c3c4tz3*conz1*tz3;
  zzcon4 = c3c4tz3*con16*tz3;
  zzcon5 = c3c4tz3*c1c5*tz3;
}
static void verify(int no_time_steps, char *class, boolean *verified) {
  double xcrref[5],xceref[5],xcrdif[5],xcedif[5], 
    epsilon, xce[5], xcr[5], dtref;
  int m;
  epsilon = 1.0e-08;
  error_norm(xce);
  compute_rhs();
  rhs_norm(xcr);
  for (m = 0; m < 5; m++) {
    xcr[m] = xcr[m] / dt;
  }
  *class = 'U';
  *verified = TRUE;
  for (m = 0; m < 5; m++) {
    xcrref[m] = 1.0;
    xceref[m] = 1.0;
  }
  if (grid_points[0] == 12 &&
      grid_points[1] == 12 &&
      grid_points[2] == 12 &&
      no_time_steps == 60) {
    *class = 'S';
    dtref = 1.0e-2;
    xcrref[0] = 1.7034283709541311e-01;
    xcrref[1] = 1.2975252070034097e-02;
    xcrref[2] = 3.2527926989486055e-02;
    xcrref[3] = 2.6436421275166801e-02;
    xcrref[4] = 1.9211784131744430e-01;
    xceref[0] = 4.9976913345811579e-04;
    xceref[1] = 4.5195666782961927e-05;
    xceref[2] = 7.3973765172921357e-05;
    xceref[3] = 7.3821238632439731e-05;
    xceref[4] = 8.9269630987491446e-04;
    } else if (grid_points[0] == 24 &&
	       grid_points[1] == 24 &&
	       grid_points[2] == 24 &&
	       no_time_steps == 200) {
      *class = 'W';
      dtref = 0.8e-3;
      xcrref[0] = 0.1125590409344e+03;
      xcrref[1] = 0.1180007595731e+02;
      xcrref[2] = 0.2710329767846e+02;
      xcrref[3] = 0.2469174937669e+02;
      xcrref[4] = 0.2638427874317e+03;
      xceref[0] = 0.4419655736008e+01;
      xceref[1] = 0.4638531260002e+00;
      xceref[2] = 0.1011551749967e+01;
      xceref[3] = 0.9235878729944e+00;
      xceref[4] = 0.1018045837718e+02;
    } else if (grid_points[0] == 64 &&
	       grid_points[1] == 64 &&
	       grid_points[2] == 64 &&
	       no_time_steps == 200) {
      *class = 'A';
      dtref = 0.8e-3;
      xcrref[0] = 1.0806346714637264e+02;
      xcrref[1] = 1.1319730901220813e+01;
      xcrref[2] = 2.5974354511582465e+01;
      xcrref[3] = 2.3665622544678910e+01;
      xcrref[4] = 2.5278963211748344e+02;
      xceref[0] = 4.2348416040525025e+00;
      xceref[1] = 4.4390282496995698e-01;
      xceref[2] = 9.6692480136345650e-01;
      xceref[3] = 8.8302063039765474e-01;
      xceref[4] = 9.7379901770829278e+00;
    } else if (grid_points[0] == 102 &&
	       grid_points[1] == 102 &&
	       grid_points[2] == 102 &&
	       no_time_steps == 200) {
      *class = 'B';
      dtref = 3.0e-4;
      xcrref[0] = 1.4233597229287254e+03;
      xcrref[1] = 9.9330522590150238e+01;
      xcrref[2] = 3.5646025644535285e+02;
      xcrref[3] = 3.2485447959084092e+02;
      xcrref[4] = 3.2707541254659363e+03;
      xceref[0] = 5.2969847140936856e+01;
      xceref[1] = 4.4632896115670668e+00;
      xceref[2] = 1.3122573342210174e+01;
      xceref[3] = 1.2006925323559144e+01;
      xceref[4] = 1.2459576151035986e+02;
    } else if (grid_points[0] == 162 &&
	       grid_points[1] == 162 &&
	       grid_points[2] == 162 &&
	       no_time_steps == 200) {
      *class = 'C';
      dtref = 1.0e-4;
      xcrref[0] = 0.62398116551764615e+04;
      xcrref[1] = 0.50793239190423964e+03;
      xcrref[2] = 0.15423530093013596e+04;
      xcrref[3] = 0.13302387929291190e+04;
      xcrref[4] = 0.11604087428436455e+05;
      xceref[0] = 0.16462008369091265e+03;
      xceref[1] = 0.11497107903824313e+02;
      xceref[2] = 0.41207446207461508e+02;
      xceref[3] = 0.37087651059694167e+02;
      xceref[4] = 0.36211053051841265e+03;
    } else {
      *verified = FALSE;
    }
  for (m = 0; m < 5; m++) {
    xcrdif[m] = fabs((xcr[m]-xcrref[m])/xcrref[m]);
    xcedif[m] = fabs((xce[m]-xceref[m])/xceref[m]);
  }
  if (*class != 'U') {
    printf(" Verification being performed for class %1c\n", *class);
    printf(" accuracy setting for epsilon = %20.13e\n", epsilon);
    if (fabs(dt-dtref) > epsilon) {
      *verified = FALSE;
      *class = 'U';
      printf(" DT does not match the reference value of %15.8e\n", dtref);
    }
  } else {
    printf(" Unknown class\n");
  }
  if (*class != 'U') {
    printf(" Comparison of RMS-norms of residual\n");
  } else {
    printf(" RMS-norms of residual\n");
  }
  for (m = 0; m < 5; m++) {
    if (*class == 'U') {
      printf("          %2d%20.13e\n", m, xcr[m]);
    } else if (xcrdif[m] > epsilon) {
      *verified = FALSE;
      printf(" FAILURE: %2d%20.13e%20.13e%20.13e\n",
	     m, xcr[m], xcrref[m], xcrdif[m]);
    } else {
      printf("          %2d%20.13e%20.13e%20.13e\n",
	     m, xcr[m], xcrref[m], xcrdif[m]);
    }
  }
  if (*class != 'U') {
    printf(" Comparison of RMS-norms of solution error\n");
  } else {
    printf(" RMS-norms of solution error\n");
  }
  for (m = 0; m < 5; m++) {
    if (*class == 'U') {
      printf("          %2d%20.13e\n", m, xce[m]);
    } else if (xcedif[m] > epsilon) {
      *verified = FALSE;
      printf(" FAILURE: %2d%20.13e%20.13e%20.13e\n",
	     m, xce[m], xceref[m], xcedif[m]);
    } else {
      printf("          %2d%20.13e%20.13e%20.13e\n",
	     m, xce[m], xceref[m], xcedif[m]);
    }
  }
  if (*class == 'U') {
    printf(" No reference values provided\n");
    printf(" No verification performed\n");
  } else if (*verified == TRUE) {
    printf(" Verification Successful\n");
  } else {
    printf(" Verification failed\n");
  }
}
static void x_solve(void) {
  lhsx();
  x_solve_cell();
  x_backsubstitute();
}
static void x_backsubstitute(void) {
  int i, j, k, m, n;
  for (i = grid_points[0]-2; i >= 0; i--) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	for (m = 0; m < BLOCK_SIZE; m++) {
	  for (n = 0; n < BLOCK_SIZE; n++) {
	    rhs[i][j][k][m] = rhs[i][j][k][m]
	      - lhs[i][j][k][CC][m][n]*rhs[i+1][j][k][n];
	  }
	}
      }
    }
  }
}
#include <omp.h> // 使用 OpenMP 需要包含这个头文件

static void x_solve_cell(void) {
  int i, j, k, isize;
  isize = grid_points[0] - 1;

  // 创建一个并行的“大扫除小组”，只创建一次，提高效率。
  // private(i, j, k) 的意思是每个组员都有自己的计数器，不会弄混。
  #pragma omp parallel private(i, j, k)
  {
    // --- 任务一：处理 i=0 的切片 ---
    // 使用 for 指令，将 j 循环的次数自动分配给小组内的成员
    #pragma omp for
    for (j = 1; j < grid_points[1] - 1; j++) {
      for (k = 1; k < grid_points[2] - 1; k++) {
        binvcrhs(lhs[0][j][k][BB],
                 lhs[0][j][k][CC],
                 rhs[0][j][k]);
      }
    }
    // omp for 结束后会有一个隐藏的“集合点”(barrier)，确保所有人都做完了再进行下一步。

    // --- 任务二：处理中间的切片 (i 从 1 到 isize-1) ---
    // i 循环必须按顺序执行，因为后一次循环依赖前一次的结果。
    for (i = 1; i < isize; i++) {
      // 对于每一个固定的 i, 内部的 j 循环可以并行。
      // 小组长说“我们现在做第 i 层”，然后大家一起上。
      #pragma omp for
      for (j = 1; j < grid_points[1] - 1; j++) {
        for (k = 1; k < grid_points[2] - 1; k++) {
          matvec_sub(lhs[i][j][k][AA],
                     rhs[i - 1][j][k], rhs[i][j][k]);
          matmul_sub(lhs[i][j][k][AA],
                     lhs[i - 1][j][k][CC],
                     lhs[i][j][k][BB]);
          binvcrhs(lhs[i][j][k][BB],
                   lhs[i][j][k][CC],
                   rhs[i][j][k]);
        }
      }
      // 每个 i 的内部任务完成后，大家又会在这里“集合”，然后再开始 i+1 的任务。
    }

    // --- 任务三：处理 i=isize 的切片 ---
    // 同样，这里的 j 循环可以并行。
    #pragma omp for
    for (j = 1; j < grid_points[1] - 1; j++) {
      for (k = 1; k < grid_points[2] - 1; k++) {
        matvec_sub(lhs[isize][j][k][AA],
                   rhs[isize - 1][j][k], rhs[isize][j][k]);
        matmul_sub(lhs[isize][j][k][AA],
                   lhs[isize - 1][j][k][CC],
                   lhs[isize][j][k][BB]);
        // 注意：原代码这里误用了变量 i, 已经修正为 isize
        binvrhs(lhs[isize][j][k][BB],
                rhs[isize][j][k]);
      }
    }
  } // “大扫除小组”工作结束
}
static void matvec_sub(double ablock[5][5], double avec[5], double bvec[5]) {
  int i;
  for (i = 0; i < 5; i++) {
    bvec[i] = bvec[i] - ablock[i][0]*avec[0]
      - ablock[i][1]*avec[1]
      - ablock[i][2]*avec[2]
      - ablock[i][3]*avec[3]
      - ablock[i][4]*avec[4];
  }
}
static void matmul_sub(double ablock[5][5], double bblock[5][5],
		       double cblock[5][5]) {
  int j;
  for (j = 0; j < 5; j++) {
    cblock[0][j] = cblock[0][j] - ablock[0][0]*bblock[0][j]
      - ablock[0][1]*bblock[1][j]
      - ablock[0][2]*bblock[2][j]
      - ablock[0][3]*bblock[3][j]
      - ablock[0][4]*bblock[4][j];
    cblock[1][j] = cblock[1][j] - ablock[1][0]*bblock[0][j]
      - ablock[1][1]*bblock[1][j]
      - ablock[1][2]*bblock[2][j]
      - ablock[1][3]*bblock[3][j]
      - ablock[1][4]*bblock[4][j];
    cblock[2][j] = cblock[2][j] - ablock[2][0]*bblock[0][j]
      - ablock[2][1]*bblock[1][j]
      - ablock[2][2]*bblock[2][j]
      - ablock[2][3]*bblock[3][j]
      - ablock[2][4]*bblock[4][j];
    cblock[3][j] = cblock[3][j] - ablock[3][0]*bblock[0][j]
      - ablock[3][1]*bblock[1][j]
      - ablock[3][2]*bblock[2][j]
      - ablock[3][3]*bblock[3][j]
      - ablock[3][4]*bblock[4][j];
    cblock[4][j] = cblock[4][j] - ablock[4][0]*bblock[0][j]
      - ablock[4][1]*bblock[1][j]
      - ablock[4][2]*bblock[2][j]
      - ablock[4][3]*bblock[3][j]
      - ablock[4][4]*bblock[4][j];
  }
}
static void binvcrhs(double lhs[5][5], double c[5][5], double r[5]) {
  double pivot, coeff;
  pivot = 1.00/lhs[0][0];
  lhs[0][1] = lhs[0][1]*pivot;
  lhs[0][2] = lhs[0][2]*pivot;
  lhs[0][3] = lhs[0][3]*pivot;
  lhs[0][4] = lhs[0][4]*pivot;
  c[0][0] = c[0][0]*pivot;
  c[0][1] = c[0][1]*pivot;
  c[0][2] = c[0][2]*pivot;
  c[0][3] = c[0][3]*pivot;
  c[0][4] = c[0][4]*pivot;
  r[0]   = r[0]  *pivot;
  coeff = lhs[1][0];
  lhs[1][1]= lhs[1][1] - coeff*lhs[0][1];
  lhs[1][2]= lhs[1][2] - coeff*lhs[0][2];
  lhs[1][3]= lhs[1][3] - coeff*lhs[0][3];
  lhs[1][4]= lhs[1][4] - coeff*lhs[0][4];
  c[1][0] = c[1][0] - coeff*c[0][0];
  c[1][1] = c[1][1] - coeff*c[0][1];
  c[1][2] = c[1][2] - coeff*c[0][2];
  c[1][3] = c[1][3] - coeff*c[0][3];
  c[1][4] = c[1][4] - coeff*c[0][4];
  r[1]   = r[1]   - coeff*r[0];
  coeff = lhs[2][0];
  lhs[2][1]= lhs[2][1] - coeff*lhs[0][1];
  lhs[2][2]= lhs[2][2] - coeff*lhs[0][2];
  lhs[2][3]= lhs[2][3] - coeff*lhs[0][3];
  lhs[2][4]= lhs[2][4] - coeff*lhs[0][4];
  c[2][0] = c[2][0] - coeff*c[0][0];
  c[2][1] = c[2][1] - coeff*c[0][1];
  c[2][2] = c[2][2] - coeff*c[0][2];
  c[2][3] = c[2][3] - coeff*c[0][3];
  c[2][4] = c[2][4] - coeff*c[0][4];
  r[2]   = r[2]   - coeff*r[0];
  coeff = lhs[3][0];
  lhs[3][1]= lhs[3][1] - coeff*lhs[0][1];
  lhs[3][2]= lhs[3][2] - coeff*lhs[0][2];
  lhs[3][3]= lhs[3][3] - coeff*lhs[0][3];
  lhs[3][4]= lhs[3][4] - coeff*lhs[0][4];
  c[3][0] = c[3][0] - coeff*c[0][0];
  c[3][1] = c[3][1] - coeff*c[0][1];
  c[3][2] = c[3][2] - coeff*c[0][2];
  c[3][3] = c[3][3] - coeff*c[0][3];
  c[3][4] = c[3][4] - coeff*c[0][4];
  r[3]   = r[3]   - coeff*r[0];
  coeff = lhs[4][0];
  lhs[4][1]= lhs[4][1] - coeff*lhs[0][1];
  lhs[4][2]= lhs[4][2] - coeff*lhs[0][2];
  lhs[4][3]= lhs[4][3] - coeff*lhs[0][3];
  lhs[4][4]= lhs[4][4] - coeff*lhs[0][4];
  c[4][0] = c[4][0] - coeff*c[0][0];
  c[4][1] = c[4][1] - coeff*c[0][1];
  c[4][2] = c[4][2] - coeff*c[0][2];
  c[4][3] = c[4][3] - coeff*c[0][3];
  c[4][4] = c[4][4] - coeff*c[0][4];
  r[4]   = r[4]   - coeff*r[0];
  pivot = 1.00/lhs[1][1];
  lhs[1][2] = lhs[1][2]*pivot;
  lhs[1][3] = lhs[1][3]*pivot;
  lhs[1][4] = lhs[1][4]*pivot;
  c[1][0] = c[1][0]*pivot;
  c[1][1] = c[1][1]*pivot;
  c[1][2] = c[1][2]*pivot;
  c[1][3] = c[1][3]*pivot;
  c[1][4] = c[1][4]*pivot;
  r[1]   = r[1]  *pivot;
  coeff = lhs[0][1];
  lhs[0][2]= lhs[0][2] - coeff*lhs[1][2];
  lhs[0][3]= lhs[0][3] - coeff*lhs[1][3];
  lhs[0][4]= lhs[0][4] - coeff*lhs[1][4];
  c[0][0] = c[0][0] - coeff*c[1][0];
  c[0][1] = c[0][1] - coeff*c[1][1];
  c[0][2] = c[0][2] - coeff*c[1][2];
  c[0][3] = c[0][3] - coeff*c[1][3];
  c[0][4] = c[0][4] - coeff*c[1][4];
  r[0]   = r[0]   - coeff*r[1];
  coeff = lhs[2][1];
  lhs[2][2]= lhs[2][2] - coeff*lhs[1][2];
  lhs[2][3]= lhs[2][3] - coeff*lhs[1][3];
  lhs[2][4]= lhs[2][4] - coeff*lhs[1][4];
  c[2][0] = c[2][0] - coeff*c[1][0];
  c[2][1] = c[2][1] - coeff*c[1][1];
  c[2][2] = c[2][2] - coeff*c[1][2];
  c[2][3] = c[2][3] - coeff*c[1][3];
  c[2][4] = c[2][4] - coeff*c[1][4];
  r[2]   = r[2]   - coeff*r[1];
  coeff = lhs[3][1];
  lhs[3][2]= lhs[3][2] - coeff*lhs[1][2];
  lhs[3][3]= lhs[3][3] - coeff*lhs[1][3];
  lhs[3][4]= lhs[3][4] - coeff*lhs[1][4];
  c[3][0] = c[3][0] - coeff*c[1][0];
  c[3][1] = c[3][1] - coeff*c[1][1];
  c[3][2] = c[3][2] - coeff*c[1][2];
  c[3][3] = c[3][3] - coeff*c[1][3];
  c[3][4] = c[3][4] - coeff*c[1][4];
  r[3]   = r[3]   - coeff*r[1];
  coeff = lhs[4][1];
  lhs[4][2]= lhs[4][2] - coeff*lhs[1][2];
  lhs[4][3]= lhs[4][3] - coeff*lhs[1][3];
  lhs[4][4]= lhs[4][4] - coeff*lhs[1][4];
  c[4][0] = c[4][0] - coeff*c[1][0];
  c[4][1] = c[4][1] - coeff*c[1][1];
  c[4][2] = c[4][2] - coeff*c[1][2];
  c[4][3] = c[4][3] - coeff*c[1][3];
  c[4][4] = c[4][4] - coeff*c[1][4];
  r[4]   = r[4]   - coeff*r[1];
  pivot = 1.00/lhs[2][2];
  lhs[2][3] = lhs[2][3]*pivot;
  lhs[2][4] = lhs[2][4]*pivot;
  c[2][0] = c[2][0]*pivot;
  c[2][1] = c[2][1]*pivot;
  c[2][2] = c[2][2]*pivot;
  c[2][3] = c[2][3]*pivot;
  c[2][4] = c[2][4]*pivot;
  r[2]   = r[2]  *pivot;
  coeff = lhs[0][2];
  lhs[0][3]= lhs[0][3] - coeff*lhs[2][3];
  lhs[0][4]= lhs[0][4] - coeff*lhs[2][4];
  c[0][0] = c[0][0] - coeff*c[2][0];
  c[0][1] = c[0][1] - coeff*c[2][1];
  c[0][2] = c[0][2] - coeff*c[2][2];
  c[0][3] = c[0][3] - coeff*c[2][3];
  c[0][4] = c[0][4] - coeff*c[2][4];
  r[0]   = r[0]   - coeff*r[2];
  coeff = lhs[1][2];
  lhs[1][3]= lhs[1][3] - coeff*lhs[2][3];
  lhs[1][4]= lhs[1][4] - coeff*lhs[2][4];
  c[1][0] = c[1][0] - coeff*c[2][0];
  c[1][1] = c[1][1] - coeff*c[2][1];
  c[1][2] = c[1][2] - coeff*c[2][2];
  c[1][3] = c[1][3] - coeff*c[2][3];
  c[1][4] = c[1][4] - coeff*c[2][4];
  r[1]   = r[1]   - coeff*r[2];
  coeff = lhs[3][2];
  lhs[3][3]= lhs[3][3] - coeff*lhs[2][3];
  lhs[3][4]= lhs[3][4] - coeff*lhs[2][4];
  c[3][0] = c[3][0] - coeff*c[2][0];
  c[3][1] = c[3][1] - coeff*c[2][1];
  c[3][2] = c[3][2] - coeff*c[2][2];
  c[3][3] = c[3][3] - coeff*c[2][3];
  c[3][4] = c[3][4] - coeff*c[2][4];
  r[3]   = r[3]   - coeff*r[2];
  coeff = lhs[4][2];
  lhs[4][3]= lhs[4][3] - coeff*lhs[2][3];
  lhs[4][4]= lhs[4][4] - coeff*lhs[2][4];
  c[4][0] = c[4][0] - coeff*c[2][0];
  c[4][1] = c[4][1] - coeff*c[2][1];
  c[4][2] = c[4][2] - coeff*c[2][2];
  c[4][3] = c[4][3] - coeff*c[2][3];
  c[4][4] = c[4][4] - coeff*c[2][4];
  r[4]   = r[4]   - coeff*r[2];
  pivot = 1.00/lhs[3][3];
  lhs[3][4] = lhs[3][4]*pivot;
  c[3][0] = c[3][0]*pivot;
  c[3][1] = c[3][1]*pivot;
  c[3][2] = c[3][2]*pivot;
  c[3][3] = c[3][3]*pivot;
  c[3][4] = c[3][4]*pivot;
  r[3]   = r[3]  *pivot;
  coeff = lhs[0][3];
  lhs[0][4]= lhs[0][4] - coeff*lhs[3][4];
  c[0][0] = c[0][0] - coeff*c[3][0];
  c[0][1] = c[0][1] - coeff*c[3][1];
  c[0][2] = c[0][2] - coeff*c[3][2];
  c[0][3] = c[0][3] - coeff*c[3][3];
  c[0][4] = c[0][4] - coeff*c[3][4];
  r[0]   = r[0]   - coeff*r[3];
  coeff = lhs[1][3];
  lhs[1][4]= lhs[1][4] - coeff*lhs[3][4];
  c[1][0] = c[1][0] - coeff*c[3][0];
  c[1][1] = c[1][1] - coeff*c[3][1];
  c[1][2] = c[1][2] - coeff*c[3][2];
  c[1][3] = c[1][3] - coeff*c[3][3];
  c[1][4] = c[1][4] - coeff*c[3][4];
  r[1]   = r[1]   - coeff*r[3];
  coeff = lhs[2][3];
  lhs[2][4]= lhs[2][4] - coeff*lhs[3][4];
  c[2][0] = c[2][0] - coeff*c[3][0];
  c[2][1] = c[2][1] - coeff*c[3][1];
  c[2][2] = c[2][2] - coeff*c[3][2];
  c[2][3] = c[2][3] - coeff*c[3][3];
  c[2][4] = c[2][4] - coeff*c[3][4];
  r[2]   = r[2]   - coeff*r[3];
  coeff = lhs[4][3];
  lhs[4][4]= lhs[4][4] - coeff*lhs[3][4];
  c[4][0] = c[4][0] - coeff*c[3][0];
  c[4][1] = c[4][1] - coeff*c[3][1];
  c[4][2] = c[4][2] - coeff*c[3][2];
  c[4][3] = c[4][3] - coeff*c[3][3];
  c[4][4] = c[4][4] - coeff*c[3][4];
  r[4]   = r[4]   - coeff*r[3];
  pivot = 1.00/lhs[4][4];
  c[4][0] = c[4][0]*pivot;
  c[4][1] = c[4][1]*pivot;
  c[4][2] = c[4][2]*pivot;
  c[4][3] = c[4][3]*pivot;
  c[4][4] = c[4][4]*pivot;
  r[4]   = r[4]  *pivot;
  coeff = lhs[0][4];
  c[0][0] = c[0][0] - coeff*c[4][0];
  c[0][1] = c[0][1] - coeff*c[4][1];
  c[0][2] = c[0][2] - coeff*c[4][2];
  c[0][3] = c[0][3] - coeff*c[4][3];
  c[0][4] = c[0][4] - coeff*c[4][4];
  r[0]   = r[0]   - coeff*r[4];
  coeff = lhs[1][4];
  c[1][0] = c[1][0] - coeff*c[4][0];
  c[1][1] = c[1][1] - coeff*c[4][1];
  c[1][2] = c[1][2] - coeff*c[4][2];
  c[1][3] = c[1][3] - coeff*c[4][3];
  c[1][4] = c[1][4] - coeff*c[4][4];
  r[1]   = r[1]   - coeff*r[4];
  coeff = lhs[2][4];
  c[2][0] = c[2][0] - coeff*c[4][0];
  c[2][1] = c[2][1] - coeff*c[4][1];
  c[2][2] = c[2][2] - coeff*c[4][2];
  c[2][3] = c[2][3] - coeff*c[4][3];
  c[2][4] = c[2][4] - coeff*c[4][4];
  r[2]   = r[2]   - coeff*r[4];
  coeff = lhs[3][4];
  c[3][0] = c[3][0] - coeff*c[4][0];
  c[3][1] = c[3][1] - coeff*c[4][1];
  c[3][2] = c[3][2] - coeff*c[4][2];
  c[3][3] = c[3][3] - coeff*c[4][3];
  c[3][4] = c[3][4] - coeff*c[4][4];
  r[3]   = r[3]   - coeff*r[4];
}
static void binvrhs( double lhs[5][5], double r[5] ) {
  double pivot, coeff;
  pivot = 1.00/lhs[0][0];
  lhs[0][1] = lhs[0][1]*pivot;
  lhs[0][2] = lhs[0][2]*pivot;
  lhs[0][3] = lhs[0][3]*pivot;
  lhs[0][4] = lhs[0][4]*pivot;
  r[0]   = r[0]  *pivot;
  coeff = lhs[1][0];
  lhs[1][1]= lhs[1][1] - coeff*lhs[0][1];
  lhs[1][2]= lhs[1][2] - coeff*lhs[0][2];
  lhs[1][3]= lhs[1][3] - coeff*lhs[0][3];
  lhs[1][4]= lhs[1][4] - coeff*lhs[0][4];
  r[1]   = r[1]   - coeff*r[0];
  coeff = lhs[2][0];
  lhs[2][1]= lhs[2][1] - coeff*lhs[0][1];
  lhs[2][2]= lhs[2][2] - coeff*lhs[0][2];
  lhs[2][3]= lhs[2][3] - coeff*lhs[0][3];
  lhs[2][4]= lhs[2][4] - coeff*lhs[0][4];
  r[2]   = r[2]   - coeff*r[0];
  coeff = lhs[3][0];
  lhs[3][1]= lhs[3][1] - coeff*lhs[0][1];
  lhs[3][2]= lhs[3][2] - coeff*lhs[0][2];
  lhs[3][3]= lhs[3][3] - coeff*lhs[0][3];
  lhs[3][4]= lhs[3][4] - coeff*lhs[0][4];
  r[3]   = r[3]   - coeff*r[0];
  coeff = lhs[4][0];
  lhs[4][1]= lhs[4][1] - coeff*lhs[0][1];
  lhs[4][2]= lhs[4][2] - coeff*lhs[0][2];
  lhs[4][3]= lhs[4][3] - coeff*lhs[0][3];
  lhs[4][4]= lhs[4][4] - coeff*lhs[0][4];
  r[4]   = r[4]   - coeff*r[0];
  pivot = 1.00/lhs[1][1];
  lhs[1][2] = lhs[1][2]*pivot;
  lhs[1][3] = lhs[1][3]*pivot;
  lhs[1][4] = lhs[1][4]*pivot;
  r[1]   = r[1]  *pivot;
  coeff = lhs[0][1];
  lhs[0][2]= lhs[0][2] - coeff*lhs[1][2];
  lhs[0][3]= lhs[0][3] - coeff*lhs[1][3];
  lhs[0][4]= lhs[0][4] - coeff*lhs[1][4];
  r[0]   = r[0]   - coeff*r[1];
  coeff = lhs[2][1];
  lhs[2][2]= lhs[2][2] - coeff*lhs[1][2];
  lhs[2][3]= lhs[2][3] - coeff*lhs[1][3];
  lhs[2][4]= lhs[2][4] - coeff*lhs[1][4];
  r[2]   = r[2]   - coeff*r[1];
  coeff = lhs[3][1];
  lhs[3][2]= lhs[3][2] - coeff*lhs[1][2];
  lhs[3][3]= lhs[3][3] - coeff*lhs[1][3];
  lhs[3][4]= lhs[3][4] - coeff*lhs[1][4];
  r[3]   = r[3]   - coeff*r[1];
  coeff = lhs[4][1];
  lhs[4][2]= lhs[4][2] - coeff*lhs[1][2];
  lhs[4][3]= lhs[4][3] - coeff*lhs[1][3];
  lhs[4][4]= lhs[4][4] - coeff*lhs[1][4];
  r[4]   = r[4]   - coeff*r[1];
  pivot = 1.00/lhs[2][2];
  lhs[2][3] = lhs[2][3]*pivot;
  lhs[2][4] = lhs[2][4]*pivot;
  r[2]   = r[2]  *pivot;
  coeff = lhs[0][2];
  lhs[0][3]= lhs[0][3] - coeff*lhs[2][3];
  lhs[0][4]= lhs[0][4] - coeff*lhs[2][4];
  r[0]   = r[0]   - coeff*r[2];
  coeff = lhs[1][2];
  lhs[1][3]= lhs[1][3] - coeff*lhs[2][3];
  lhs[1][4]= lhs[1][4] - coeff*lhs[2][4];
  r[1]   = r[1]   - coeff*r[2];
  coeff = lhs[3][2];
  lhs[3][3]= lhs[3][3] - coeff*lhs[2][3];
  lhs[3][4]= lhs[3][4] - coeff*lhs[2][4];
  r[3]   = r[3]   - coeff*r[2];
  coeff = lhs[4][2];
  lhs[4][3]= lhs[4][3] - coeff*lhs[2][3];
  lhs[4][4]= lhs[4][4] - coeff*lhs[2][4];
  r[4]   = r[4]   - coeff*r[2];
  pivot = 1.00/lhs[3][3];
  lhs[3][4] = lhs[3][4]*pivot;
  r[3]   = r[3]  *pivot;
  coeff = lhs[0][3];
  lhs[0][4]= lhs[0][4] - coeff*lhs[3][4];
  r[0]   = r[0]   - coeff*r[3];
  coeff = lhs[1][3];
  lhs[1][4]= lhs[1][4] - coeff*lhs[3][4];
  r[1]   = r[1]   - coeff*r[3];
  coeff = lhs[2][3];
  lhs[2][4]= lhs[2][4] - coeff*lhs[3][4];
  r[2]   = r[2]   - coeff*r[3];
  coeff = lhs[4][3];
  lhs[4][4]= lhs[4][4] - coeff*lhs[3][4];
  r[4]   = r[4]   - coeff*r[3];
  pivot = 1.00/lhs[4][4];
  r[4]   = r[4]  *pivot;
  coeff = lhs[0][4];
  r[0]   = r[0]   - coeff*r[4];
  coeff = lhs[1][4];
  r[1]   = r[1]   - coeff*r[4];
  coeff = lhs[2][4];
  r[2]   = r[2]   - coeff*r[4];
  coeff = lhs[3][4];
  r[3]   = r[3]   - coeff*r[4];
}
static void y_solve(void) {
  lhsy();
  y_solve_cell();
  y_backsubstitute();
}
static void y_backsubstitute(void) {
  int i, j, k, m, n;
  for (j = grid_points[1]-2; j >= 0; j--) {
    for (i = 1; i < grid_points[0]-1; i++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	for (m = 0; m < BLOCK_SIZE; m++) {
	  for (n = 0; n < BLOCK_SIZE; n++) {
	    rhs[i][j][k][m] = rhs[i][j][k][m] 
	      - lhs[i][j][k][CC][m][n]*rhs[i][j+1][k][n];
	  }
	}
      }
    }
  }
}
#include <omp.h>

static void y_solve_cell(void) {
  int i, j, k, jsize;

  // 使用一个 parallel 区域包裹整个函数，避免重复创建/销毁线程
  // i, j, k 作为循环变量，必须是线程私有的（private）
  // jsize 也在并行区域内定义，因此也是线程私有的
  #pragma omp parallel private(i, j, k, jsize)
  {
    // ---------------------------------------------------------------------
    // 第一部分: 初始化边界条件 (j=0)
    // ---------------------------------------------------------------------
    // i 和 k 循环的迭代之间没有依赖关系，可以安全地并行化。
    // 使用 collapse(2) 将两个循环合并，以获得更好的负载均衡。
    #pragma omp for collapse(2)
    for (i = 1; i < grid_points[0]-1; i++) {
      for (k = 1; k < grid_points[2]-1; k++) {
        binvcrhs( lhs[i][0][k][BB],
                  lhs[i][0][k][CC],
                  rhs[i][0][k] );
      }
    }
    // 在 #pragma omp for 之后有一个隐式屏障。
    // 所有线程必须完成这里的计算，才能继续执行下面的代码。
    // 这保证了 j=0 的结果在 j=1 开始计算前已经准备好。

    // ---------------------------------------------------------------------
    // 第二部分: 前向替换的主循环
    // ---------------------------------------------------------------------
    // j 循环有前后依赖关系 (j 依赖 j-1)，因此必须串行执行。
    // 但是在每次 j 迭代的内部，i 和 k 的计算是独立的，可以并行化。
    jsize = grid_points[1]-1;
    for (j = 1; j < jsize; j++) {
      // 在串行 j 循环内部，将 i-k 循环的工作分配给所有线程。
      #pragma omp for collapse(2)
      for (i = 1; i < grid_points[0]-1; i++) {
        for (k = 1; k < grid_points[2]-1; k++) {
          matvec_sub(lhs[i][j][k][AA],
                     rhs[i][j-1][k], rhs[i][j][k]);
          matmul_sub(lhs[i][j][k][AA],
                     lhs[i][j-1][k][CC],
                     lhs[i][j][k][BB]);
          binvcrhs( lhs[i][j][k][BB],
                    lhs[i][j][k][CC],
                    rhs[i][j][k] );
        }
      }
      // 同样，在每次 #pragma omp for 之后都有一个隐式屏障。
      // 这确保了对于给定的 j，所有 (i, k) 的计算都已完成，
      // 之后主线程（逻辑上）才会进入 j+1 的迭代。
    }

    // ---------------------------------------------------------------------
    // 第三部分: 回代求解 (j=jsize)
    // ---------------------------------------------------------------------
    // 主循环完成后，所有线程已同步，j=jsize-1 的结果已可用。
    // 这里的 i-k 循环同样没有内部依赖，可以并行化。
    #pragma omp for collapse(2)
    for (i = 1; i < grid_points[0]-1; i++) {
      for (k = 1; k < grid_points[2]-1; k++) {
        matvec_sub(lhs[i][jsize][k][AA],
                   rhs[i][jsize-1][k], rhs[i][jsize][k]);
        matmul_sub(lhs[i][jsize][k][AA],
                   lhs[i][jsize-1][k][CC],
                   lhs[i][jsize][k][BB]);
        binvrhs( lhs[i][jsize][k][BB],
                 rhs[i][jsize][k] );
      }
    }
    // 并行区域末尾有隐式屏障，所有线程在此同步后退出。
  } // end of #pragma omp parallel
}
static void z_solve(void) {
  lhsz();
  z_solve_cell();
  z_backsubstitute();
}
static void z_backsubstitute(void) {
  int i, j, k, m, n;
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = grid_points[2]-2; k >= 0; k--) {
	for (m = 0; m < BLOCK_SIZE; m++) {
	  for (n = 0; n < BLOCK_SIZE; n++) {
	    rhs[i][j][k][m] = rhs[i][j][k][m] 
	      - lhs[i][j][k][CC][m][n]*rhs[i][j][k+1][n];
	  }
	}
      }
    }
  }
}
#include <omp.h>


static void z_solve_cell(void) {
  int i, j, k, ksize;
  ksize = grid_points[2] - 1;

  // 使用单一的并行区域来减少线程创建/销毁的开销
  #pragma omp parallel private(i, j, k)
  {
    //---------------------------------------------------------------------
    // 第一个循环: 对 k=0 平面进行计算
    // i,j 循环是独立的，可以并行化
    // 使用 collapse(2) 将两层循环合并为一个大的迭代空间，可以改善负载均衡
    //---------------------------------------------------------------------
    #pragma omp for collapse(2)
    for (i = 1; i < grid_points[0] - 1; i++) {
      for (j = 1; j < grid_points[1] - 1; j++) {
        binvcrhs(lhs[i][j][0][BB], 
                 lhs[i][j][0][CC], 
                 rhs[i][j][0]);
      }
    }
    // omp for 结束时有一个隐式屏障，确保所有线程都完成了 k=0 的计算

    //---------------------------------------------------------------------
    // 第二个循环: 前向替换过程
    // k 循环存在依赖性 (k 依赖 k-1)，必须串行执行
    // 对于每个 k，内部的 i,j 循环是独立的，可以并行化
    //---------------------------------------------------------------------
    for (k = 1; k < ksize; k++) {
      // 在串行的 k 循环内部，并行化 i 和 j 循环
      #pragma omp for collapse(2)
      for (i = 1; i < grid_points[0] - 1; i++) {
        for (j = 1; j < grid_points[1] - 1; j++) {
          matvec_sub(lhs[i][j][k][AA], 
                     rhs[i][j][k - 1], rhs[i][j][k]);
          matmul_sub(lhs[i][j][k][AA], 
                     lhs[i][j][k - 1][CC], 
                     lhs[i][j][k][BB]);
          binvcrhs(lhs[i][j][k][BB], 
                   lhs[i][j][k][CC], 
                   rhs[i][j][k]);
        }
      }
      // omp for 结束时有一个隐式屏障
      // 这确保了在进入下一次 k+1 迭代之前，当前 k 的所有计算都已完成
    }

    //---------------------------------------------------------------------
    // 第三个循环: 对 k=ksize 平面进行计算
    // i,j 循环是独立的，可以并行化
    //---------------------------------------------------------------------
    #pragma omp for collapse(2)
    for (i = 1; i < grid_points[0] - 1; i++) {
      for (j = 1; j < grid_points[1] - 1; j++) {
        matvec_sub(lhs[i][j][ksize][AA], 
                   rhs[i][j][ksize - 1], rhs[i][j][ksize]);
        matmul_sub(lhs[i][j][ksize][AA], 
                   lhs[i][j][ksize - 1][CC],
                   lhs[i][j][ksize][BB]);
        binvrhs(lhs[i][j][ksize][BB], 
                rhs[i][j][ksize]);
      }
    }
    // omp for 结束时有隐式屏障，之后整个并行区域结束
  } // pragma omp parallel
}