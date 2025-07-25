#include "npb-C.h"
#include "header.h"
#include <omp.h> 
static void add();
static void adi();
static void error_norm(double rms[5]);
static void rhs_norm(double rms[5]);
static void exact_rhs();
static void exact_solution(double xi,double eta,double zeta,double dtemp[5]);
static void initialize();
static void lhsinit();
static void lhsx();
static void lhsy();
static void lhsz();
static void compute_rhs();
static void set_constants();
static void verify(int no_time_steps,char *class,boolean *verified);
static void x_solve();
static void x_backsubstitute();
static void x_solve_cell();
static void matvec_sub(double ablock[5][5],double avec[5],double bvec[5]);
static void matmul_sub(double ablock[5][5],double bblock[5][5],double cblock[5][5]);
static void binvcrhs(double lhs[5][5],double c[5][5],double r[5]);
static void binvrhs(double lhs[5][5],double r[5]);
static void y_solve();
static void y_backsubstitute();
static void y_solve_cell();
static void z_solve();
static void z_backsubstitute();
static void z_solve_cell();

int main(int argc,char **argv)
{
  int niter;
  int step;
  int n3;
  int nthreads = 1;
  double navg;
  double mflops;
  double tmax;
  boolean verified;
  char class;
  FILE *fp;
  printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version - BT Benchmark\n\n");
  fp = fopen("inputbt.data","r");
  if (fp != ((void *)0)) {
    printf(" Reading from input file inputbt.data");
    fscanf(fp,"%d",&niter);
    while(fgetc(fp) != '\n')
      ;
    fscanf(fp,"%lg",&dt);
    while(fgetc(fp) != '\n')
      ;
    fscanf(fp,"%d%d%d",&grid_points[0],&grid_points[1],&grid_points[2]);
    fclose(fp);
  }
   else {
    printf(" No input file inputbt.data. Using compiled defaults\n");
    niter = 200;
    dt = 0.0008;
    grid_points[0] = 24;
    grid_points[1] = 24;
    grid_points[2] = 24;
  }
  printf(" Size: %3dx%3dx%3d\n",grid_points[0],grid_points[1],grid_points[2]);
  printf(" Iterations: %3d   dt: %10.6f\n",niter,dt);
  if (grid_points[0] > 24 || grid_points[1] > 24 || grid_points[2] > 24) {
    printf(" %dx%dx%d\n",grid_points[0],grid_points[1],grid_points[2]);
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
    for (step = 1; step <= niter; step += 1) {
      if (step % 20 == 0 || step == 1) {
        printf(" Time step %4d\n",step);
      }
      adi();
    }
  }
{
    #if defined(_OPENMP)
    #endif 
  }
  timer_stop(1);
  tmax = timer_read(1);
{
    verify(niter,&class,&verified);
  }
  n3 = grid_points[0] * grid_points[1] * grid_points[2];
  navg = (grid_points[0] + grid_points[1] + grid_points[2]) / 3.0;
  if (tmax != 0.0) {
    mflops = 1.0e-6 * ((double )niter) * (3478.8 * ((double )n3) - 17655.7 * (navg * navg) + 28023.7 * navg) / tmax;
  }
   else {
    mflops = 0.0;
  }
  c_print_results("BT",class,grid_points[0],grid_points[1],grid_points[2],niter,nthreads,tmax,mflops,"          floating point",verified,"3.0 structured","25 Jun 2025","gcc","gcc","-lm","-I../common","-O3 -fopenmp -fopt-info-vec-missed=vec-miss...","-fopenmp -lm","(none)");
}

static void add()
{
  int i;
  int j;
  int k;
  int m;
  
#pragma omp parallel for private (i,j,k,m)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (j,k,m)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      
#pragma omp parallel for private (k,m)
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        
#pragma omp parallel for private (m)
        for (m = 0; m <= 4; m += 1) {
          u[i][j][k][m] = u[i][j][k][m] + rhs[i][j][k][m];
        }
      }
    }
  }
}

static void adi()
{
  compute_rhs();
  x_solve();
  y_solve();
  z_solve();
  add();
}

static void error_norm(double rms[5])
{
  int i;
  int j;
  int k;
  int m;
  int d;
  double xi;
  double eta;
  double zeta;
  double u_exact[5];
  double add;
  
#pragma omp parallel for private (m)
  for (m = 0; m <= 4; m += 1) {
    rms[m] = 0.0;
  }
  for (i = 0; i <= grid_points[0] - 1; i += 1) {
    xi = ((double )i) * dnxm1;
    for (j = 0; j <= grid_points[1] - 1; j += 1) {
      eta = ((double )j) * dnym1;
      for (k = 0; k <= grid_points[2] - 1; k += 1) {
        zeta = ((double )k) * dnzm1;
        exact_solution(xi,eta,zeta,u_exact);
        
#pragma omp parallel for private (add,m)
        for (m = 0; m <= 4; m += 1) {
          add = u[i][j][k][m] - u_exact[m];
          rms[m] = rms[m] + add * add;
        }
      }
    }
  }
  for (m = 0; m <= 4; m += 1) {
    for (d = 0; d <= 2; d += 1) {
      rms[m] = rms[m] / ((double )(grid_points[d] - 2));
    }
    rms[m] = sqrt(rms[m]);
  }
}

static void rhs_norm(double rms[5])
{
  int i;
  int j;
  int k;
  int d;
  int m;
  double add;
  
#pragma omp parallel for private (m)
  for (m = 0; m <= 4; m += 1) {
    rms[m] = 0.0;
  }
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        
#pragma omp parallel for private (add,m)
        for (m = 0; m <= 4; m += 1) {
          add = rhs[i][j][k][m];
          rms[m] = rms[m] + add * add;
        }
      }
    }
  }
  for (m = 0; m <= 4; m += 1) {
    for (d = 0; d <= 2; d += 1) {
      rms[m] = rms[m] / ((double )(grid_points[d] - 2));
    }
    rms[m] = sqrt(rms[m]);
  }
}

static void exact_rhs()
{
{
    double dtemp[5];
    double xi;
    double eta;
    double zeta;
    double dtpp;
    int m;
    int i;
    int j;
    int k;
    int ip1;
    int im1;
    int jp1;
    int jm1;
    int km1;
    int kp1;
    
#pragma omp parallel for private (m,i,j,k)
    for (i = 0; i <= grid_points[0] - 1; i += 1) {
      
#pragma omp parallel for private (m,j,k)
      for (j = 0; j <= grid_points[1] - 1; j += 1) {
        
#pragma omp parallel for private (m,k)
        for (k = 0; k <= grid_points[2] - 1; k += 1) {
          
#pragma omp parallel for private (m)
          for (m = 0; m <= 4; m += 1) {
            forcing[i][j][k][m] = 0.0;
          }
        }
      }
    }
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      eta = ((double )j) * dnym1;
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        zeta = ((double )k) * dnzm1;
        for (i = 0; i <= grid_points[0] - 1; i += 1) {
          xi = ((double )i) * dnxm1;
          exact_solution(xi,eta,zeta,dtemp);
          
#pragma omp parallel for private (m)
          for (m = 0; m <= 4; m += 1) {
            ue[i][m] = dtemp[m];
          }
          dtpp = 1.0 / dtemp[0];
          
#pragma omp parallel for private (m) firstprivate (dtpp)
          for (m = 1; m <= 4; m += 1) {
            buf[i][m] = dtpp * dtemp[m];
          }
          cuf[i] = buf[i][1] * buf[i][1];
          buf[i][0] = cuf[i] + buf[i][2] * buf[i][2] + buf[i][3] * buf[i][3];
          q[i] = 0.5 * (buf[i][1] * ue[i][1] + buf[i][2] * ue[i][2] + buf[i][3] * ue[i][3]);
        }
        for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
          im1 = i - 1;
          ip1 = i + 1;
          forcing[i][j][k][0] = forcing[i][j][k][0] - tx2 * (ue[ip1][1] - ue[im1][1]) + dx1tx1 * (ue[ip1][0] - 2.0 * ue[i][0] + ue[im1][0]);
          forcing[i][j][k][1] = forcing[i][j][k][1] - tx2 * (ue[ip1][1] * buf[ip1][1] + c2 * (ue[ip1][4] - q[ip1]) - (ue[im1][1] * buf[im1][1] + c2 * (ue[im1][4] - q[im1]))) + xxcon1 * (buf[ip1][1] - 2.0 * buf[i][1] + buf[im1][1]) + dx2tx1 * (ue[ip1][1] - 2.0 * ue[i][1] + ue[im1][1]);
          forcing[i][j][k][2] = forcing[i][j][k][2] - tx2 * (ue[ip1][2] * buf[ip1][1] - ue[im1][2] * buf[im1][1]) + xxcon2 * (buf[ip1][2] - 2.0 * buf[i][2] + buf[im1][2]) + dx3tx1 * (ue[ip1][2] - 2.0 * ue[i][2] + ue[im1][2]);
          forcing[i][j][k][3] = forcing[i][j][k][3] - tx2 * (ue[ip1][3] * buf[ip1][1] - ue[im1][3] * buf[im1][1]) + xxcon2 * (buf[ip1][3] - 2.0 * buf[i][3] + buf[im1][3]) + dx4tx1 * (ue[ip1][3] - 2.0 * ue[i][3] + ue[im1][3]);
          forcing[i][j][k][4] = forcing[i][j][k][4] - tx2 * (buf[ip1][1] * (c1 * ue[ip1][4] - c2 * q[ip1]) - buf[im1][1] * (c1 * ue[im1][4] - c2 * q[im1])) + 0.5 * xxcon3 * (buf[ip1][0] - 2.0 * buf[i][0] + buf[im1][0]) + xxcon4 * (cuf[ip1] - 2.0 * cuf[i] + cuf[im1]) + xxcon5 * (buf[ip1][4] - 2.0 * buf[i][4] + buf[im1][4]) + dx5tx1 * (ue[ip1][4] - 2.0 * ue[i][4] + ue[im1][4]);
        }
        
#pragma omp parallel for private (i,m)
        for (m = 0; m <= 4; m += 1) {
          i = 1;
          forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (5.0 * ue[i][m] - 4.0 * ue[i + 1][m] + ue[i + 2][m]);
          i = 2;
          forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (- 4.0 * ue[i - 1][m] + 6.0 * ue[i][m] - 4.0 * ue[i + 1][m] + ue[i + 2][m]);
        }
        
#pragma omp parallel for private (m,i)
        for (m = 0; m <= 4; m += 1) {
          for (i = 1 * 3; i <= grid_points[0] - 3 - 1; i += 1) {
            forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[i - 2][m] - 4.0 * ue[i - 1][m] + 6.0 * ue[i][m] - 4.0 * ue[i + 1][m] + ue[i + 2][m]);
          }
        }
        for (m = 0; m <= 4; m += 1) {
          i = grid_points[0] - 3;
          forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[i - 2][m] - 4.0 * ue[i - 1][m] + 6.0 * ue[i][m] - 4.0 * ue[i + 1][m]);
          i = grid_points[0] - 2;
          forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[i - 2][m] - 4.0 * ue[i - 1][m] + 5.0 * ue[i][m]);
        }
      }
    }
    for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
      xi = ((double )i) * dnxm1;
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        zeta = ((double )k) * dnzm1;
        for (j = 0; j <= grid_points[1] - 1; j += 1) {
          eta = ((double )j) * dnym1;
          exact_solution(xi,eta,zeta,dtemp);
          
#pragma omp parallel for private (m)
          for (m = 0; m <= 4; m += 1) {
            ue[j][m] = dtemp[m];
          }
          dtpp = 1.0 / dtemp[0];
          
#pragma omp parallel for private (m) firstprivate (dtpp)
          for (m = 1; m <= 4; m += 1) {
            buf[j][m] = dtpp * dtemp[m];
          }
          cuf[j] = buf[j][2] * buf[j][2];
          buf[j][0] = cuf[j] + buf[j][1] * buf[j][1] + buf[j][3] * buf[j][3];
          q[j] = 0.5 * (buf[j][1] * ue[j][1] + buf[j][2] * ue[j][2] + buf[j][3] * ue[j][3]);
        }
        for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
          jm1 = j - 1;
          jp1 = j + 1;
          forcing[i][j][k][0] = forcing[i][j][k][0] - ty2 * (ue[jp1][2] - ue[jm1][2]) + dy1ty1 * (ue[jp1][0] - 2.0 * ue[j][0] + ue[jm1][0]);
          forcing[i][j][k][1] = forcing[i][j][k][1] - ty2 * (ue[jp1][1] * buf[jp1][2] - ue[jm1][1] * buf[jm1][2]) + yycon2 * (buf[jp1][1] - 2.0 * buf[j][1] + buf[jm1][1]) + dy2ty1 * (ue[jp1][1] - 2.0 * ue[j][1] + ue[jm1][1]);
          forcing[i][j][k][2] = forcing[i][j][k][2] - ty2 * (ue[jp1][2] * buf[jp1][2] + c2 * (ue[jp1][4] - q[jp1]) - (ue[jm1][2] * buf[jm1][2] + c2 * (ue[jm1][4] - q[jm1]))) + yycon1 * (buf[jp1][2] - 2.0 * buf[j][2] + buf[jm1][2]) + dy3ty1 * (ue[jp1][2] - 2.0 * ue[j][2] + ue[jm1][2]);
          forcing[i][j][k][3] = forcing[i][j][k][3] - ty2 * (ue[jp1][3] * buf[jp1][2] - ue[jm1][3] * buf[jm1][2]) + yycon2 * (buf[jp1][3] - 2.0 * buf[j][3] + buf[jm1][3]) + dy4ty1 * (ue[jp1][3] - 2.0 * ue[j][3] + ue[jm1][3]);
          forcing[i][j][k][4] = forcing[i][j][k][4] - ty2 * (buf[jp1][2] * (c1 * ue[jp1][4] - c2 * q[jp1]) - buf[jm1][2] * (c1 * ue[jm1][4] - c2 * q[jm1])) + 0.5 * yycon3 * (buf[jp1][0] - 2.0 * buf[j][0] + buf[jm1][0]) + yycon4 * (cuf[jp1] - 2.0 * cuf[j] + cuf[jm1]) + yycon5 * (buf[jp1][4] - 2.0 * buf[j][4] + buf[jm1][4]) + dy5ty1 * (ue[jp1][4] - 2.0 * ue[j][4] + ue[jm1][4]);
        }
        
#pragma omp parallel for private (j,m)
        for (m = 0; m <= 4; m += 1) {
          j = 1;
          forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (5.0 * ue[j][m] - 4.0 * ue[j + 1][m] + ue[j + 2][m]);
          j = 2;
          forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (- 4.0 * ue[j - 1][m] + 6.0 * ue[j][m] - 4.0 * ue[j + 1][m] + ue[j + 2][m]);
        }
        
#pragma omp parallel for private (m,j)
        for (m = 0; m <= 4; m += 1) {
          for (j = 1 * 3; j <= grid_points[1] - 3 - 1; j += 1) {
            forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[j - 2][m] - 4.0 * ue[j - 1][m] + 6.0 * ue[j][m] - 4.0 * ue[j + 1][m] + ue[j + 2][m]);
          }
        }
        for (m = 0; m <= 4; m += 1) {
          j = grid_points[1] - 3;
          forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[j - 2][m] - 4.0 * ue[j - 1][m] + 6.0 * ue[j][m] - 4.0 * ue[j + 1][m]);
          j = grid_points[1] - 2;
          forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[j - 2][m] - 4.0 * ue[j - 1][m] + 5.0 * ue[j][m]);
        }
      }
    }
    for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
      xi = ((double )i) * dnxm1;
      for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
        eta = ((double )j) * dnym1;
        for (k = 0; k <= grid_points[2] - 1; k += 1) {
          zeta = ((double )k) * dnzm1;
          exact_solution(xi,eta,zeta,dtemp);
          
#pragma omp parallel for private (m)
          for (m = 0; m <= 4; m += 1) {
            ue[k][m] = dtemp[m];
          }
          dtpp = 1.0 / dtemp[0];
          
#pragma omp parallel for private (m) firstprivate (dtpp)
          for (m = 1; m <= 4; m += 1) {
            buf[k][m] = dtpp * dtemp[m];
          }
          cuf[k] = buf[k][3] * buf[k][3];
          buf[k][0] = cuf[k] + buf[k][1] * buf[k][1] + buf[k][2] * buf[k][2];
          q[k] = 0.5 * (buf[k][1] * ue[k][1] + buf[k][2] * ue[k][2] + buf[k][3] * ue[k][3]);
        }
        for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
          km1 = k - 1;
          kp1 = k + 1;
          forcing[i][j][k][0] = forcing[i][j][k][0] - tz2 * (ue[kp1][3] - ue[km1][3]) + dz1tz1 * (ue[kp1][0] - 2.0 * ue[k][0] + ue[km1][0]);
          forcing[i][j][k][1] = forcing[i][j][k][1] - tz2 * (ue[kp1][1] * buf[kp1][3] - ue[km1][1] * buf[km1][3]) + zzcon2 * (buf[kp1][1] - 2.0 * buf[k][1] + buf[km1][1]) + dz2tz1 * (ue[kp1][1] - 2.0 * ue[k][1] + ue[km1][1]);
          forcing[i][j][k][2] = forcing[i][j][k][2] - tz2 * (ue[kp1][2] * buf[kp1][3] - ue[km1][2] * buf[km1][3]) + zzcon2 * (buf[kp1][2] - 2.0 * buf[k][2] + buf[km1][2]) + dz3tz1 * (ue[kp1][2] - 2.0 * ue[k][2] + ue[km1][2]);
          forcing[i][j][k][3] = forcing[i][j][k][3] - tz2 * (ue[kp1][3] * buf[kp1][3] + c2 * (ue[kp1][4] - q[kp1]) - (ue[km1][3] * buf[km1][3] + c2 * (ue[km1][4] - q[km1]))) + zzcon1 * (buf[kp1][3] - 2.0 * buf[k][3] + buf[km1][3]) + dz4tz1 * (ue[kp1][3] - 2.0 * ue[k][3] + ue[km1][3]);
          forcing[i][j][k][4] = forcing[i][j][k][4] - tz2 * (buf[kp1][3] * (c1 * ue[kp1][4] - c2 * q[kp1]) - buf[km1][3] * (c1 * ue[km1][4] - c2 * q[km1])) + 0.5 * zzcon3 * (buf[kp1][0] - 2.0 * buf[k][0] + buf[km1][0]) + zzcon4 * (cuf[kp1] - 2.0 * cuf[k] + cuf[km1]) + zzcon5 * (buf[kp1][4] - 2.0 * buf[k][4] + buf[km1][4]) + dz5tz1 * (ue[kp1][4] - 2.0 * ue[k][4] + ue[km1][4]);
        }
        
#pragma omp parallel for private (k,m)
        for (m = 0; m <= 4; m += 1) {
          k = 1;
          forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (5.0 * ue[k][m] - 4.0 * ue[k + 1][m] + ue[k + 2][m]);
          k = 2;
          forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (- 4.0 * ue[k - 1][m] + 6.0 * ue[k][m] - 4.0 * ue[k + 1][m] + ue[k + 2][m]);
        }
        
#pragma omp parallel for private (m,k)
        for (m = 0; m <= 4; m += 1) {
          for (k = 1 * 3; k <= grid_points[2] - 3 - 1; k += 1) {
            forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[k - 2][m] - 4.0 * ue[k - 1][m] + 6.0 * ue[k][m] - 4.0 * ue[k + 1][m] + ue[k + 2][m]);
          }
        }
        for (m = 0; m <= 4; m += 1) {
          k = grid_points[2] - 3;
          forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[k - 2][m] - 4.0 * ue[k - 1][m] + 6.0 * ue[k][m] - 4.0 * ue[k + 1][m]);
          k = grid_points[2] - 2;
          forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[k - 2][m] - 4.0 * ue[k - 1][m] + 5.0 * ue[k][m]);
        }
      }
    }
    
#pragma omp parallel for private (m,i,j,k)
    for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
      
#pragma omp parallel for private (m,j,k)
      for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
        
#pragma omp parallel for private (m,k)
        for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
          
#pragma omp parallel for private (m)
          for (m = 0; m <= 4; m += 1) {
            forcing[i][j][k][m] = - 1.0 * forcing[i][j][k][m];
          }
        }
      }
    }
  }
}

static void exact_solution(double xi,double eta,double zeta,double dtemp[5])
{
  int m;
  
#pragma omp parallel for private (m) firstprivate (xi,eta,zeta)
  for (m = 0; m <= 4; m += 1) {
    dtemp[m] = ce[m][0] + xi * (ce[m][1] + xi * (ce[m][4] + xi * (ce[m][7] + xi * ce[m][10]))) + eta * (ce[m][2] + eta * (ce[m][5] + eta * (ce[m][8] + eta * ce[m][11]))) + zeta * (ce[m][3] + zeta * (ce[m][6] + zeta * (ce[m][9] + zeta * ce[m][12])));
  }
}

static void initialize()
{
{
    int i;
    int j;
    int k;
    int m;
    int ix;
    int iy;
    int iz;
    double xi;
    double eta;
    double zeta;
    double Pface[2][3][5];
    double Pxi;
    double Peta;
    double Pzeta;
    double temp[5];
    
#pragma omp parallel for private (i,j,k,m)
    for (i = 0; i <= 23; i += 1) {
      
#pragma omp parallel for private (j,k,m)
      for (j = 0; j <= 23; j += 1) {
        
#pragma omp parallel for private (k,m)
        for (k = 0; k <= 23; k += 1) {
          
#pragma omp parallel for private (m)
          for (m = 0; m <= 4; m += 1) {
            u[i][j][k][m] = 1.0;
          }
        }
      }
    }
    for (i = 0; i <= grid_points[0] - 1; i += 1) {
      xi = ((double )i) * dnxm1;
      for (j = 0; j <= grid_points[1] - 1; j += 1) {
        eta = ((double )j) * dnym1;
        for (k = 0; k <= grid_points[2] - 1; k += 1) {
          zeta = ((double )k) * dnzm1;
          for (ix = 0; ix <= 1; ix += 1) {
            exact_solution((double )ix,eta,zeta,&Pface[ix][0][0]);
          }
          for (iy = 0; iy <= 1; iy += 1) {
            exact_solution(xi,(double )iy,zeta,&Pface[iy][1][0]);
          }
          for (iz = 0; iz <= 1; iz += 1) {
            exact_solution(xi,eta,(double )iz,&Pface[iz][2][0]);
          }
          
#pragma omp parallel for private (Pxi,Peta,Pzeta,m) firstprivate (xi,eta,zeta)
          for (m = 0; m <= 4; m += 1) {
            Pxi = xi * Pface[1][0][m] + (1.0 - xi) * Pface[0][0][m];
            Peta = eta * Pface[1][1][m] + (1.0 - eta) * Pface[0][1][m];
            Pzeta = zeta * Pface[1][2][m] + (1.0 - zeta) * Pface[0][2][m];
            u[i][j][k][m] = Pxi + Peta + Pzeta - Pxi * Peta - Pxi * Pzeta - Peta * Pzeta + Pxi * Peta * Pzeta;
          }
        }
      }
    }
    i = 0;
    xi = 0.0;
    for (j = 0; j <= grid_points[1] - 1; j += 1) {
      eta = ((double )j) * dnym1;
      for (k = 0; k <= grid_points[2] - 1; k += 1) {
        zeta = ((double )k) * dnzm1;
        exact_solution(xi,eta,zeta,temp);
        
#pragma omp parallel for private (m) firstprivate (i)
        for (m = 0; m <= 4; m += 1) {
          u[i][j][k][m] = temp[m];
        }
      }
    }
    i = grid_points[0] - 1;
    xi = 1.0;
    for (j = 0; j <= grid_points[1] - 1; j += 1) {
      eta = ((double )j) * dnym1;
      for (k = 0; k <= grid_points[2] - 1; k += 1) {
        zeta = ((double )k) * dnzm1;
        exact_solution(xi,eta,zeta,temp);
        
#pragma omp parallel for private (m) firstprivate (i)
        for (m = 0; m <= 4; m += 1) {
          u[i][j][k][m] = temp[m];
        }
      }
    }
    j = 0;
    eta = 0.0;
    for (i = 0; i <= grid_points[0] - 1; i += 1) {
      xi = ((double )i) * dnxm1;
      for (k = 0; k <= grid_points[2] - 1; k += 1) {
        zeta = ((double )k) * dnzm1;
        exact_solution(xi,eta,zeta,temp);
        
#pragma omp parallel for private (m) firstprivate (j)
        for (m = 0; m <= 4; m += 1) {
          u[i][j][k][m] = temp[m];
        }
      }
    }
    j = grid_points[1] - 1;
    eta = 1.0;
    for (i = 0; i <= grid_points[0] - 1; i += 1) {
      xi = ((double )i) * dnxm1;
      for (k = 0; k <= grid_points[2] - 1; k += 1) {
        zeta = ((double )k) * dnzm1;
        exact_solution(xi,eta,zeta,temp);
        
#pragma omp parallel for private (m) firstprivate (j)
        for (m = 0; m <= 4; m += 1) {
          u[i][j][k][m] = temp[m];
        }
      }
    }
    k = 0;
    zeta = 0.0;
    for (i = 0; i <= grid_points[0] - 1; i += 1) {
      xi = ((double )i) * dnxm1;
      for (j = 0; j <= grid_points[1] - 1; j += 1) {
        eta = ((double )j) * dnym1;
        exact_solution(xi,eta,zeta,temp);
        
#pragma omp parallel for private (m) firstprivate (k)
        for (m = 0; m <= 4; m += 1) {
          u[i][j][k][m] = temp[m];
        }
      }
    }
    k = grid_points[2] - 1;
    zeta = 1.0;
    for (i = 0; i <= grid_points[0] - 1; i += 1) {
      xi = ((double )i) * dnxm1;
      for (j = 0; j <= grid_points[1] - 1; j += 1) {
        eta = ((double )j) * dnym1;
        exact_solution(xi,eta,zeta,temp);
        
#pragma omp parallel for private (m) firstprivate (k)
        for (m = 0; m <= 4; m += 1) {
          u[i][j][k][m] = temp[m];
        }
      }
    }
  }
}

static void lhsinit()
{
{
    int i;
    int j;
    int k;
    int m;
    int n;
    
#pragma omp parallel for private (i,j,k,m,n)
    for (i = 0; i <= grid_points[0] - 1; i += 1) {
      
#pragma omp parallel for private (j,k,m,n)
      for (j = 0; j <= grid_points[1] - 1; j += 1) {
        
#pragma omp parallel for private (k,m,n)
        for (k = 0; k <= grid_points[2] - 1; k += 1) {
          
#pragma omp parallel for private (m,n)
          for (m = 0; m <= 4; m += 1) {
            
#pragma omp parallel for private (n)
            for (n = 0; n <= 4; n += 1) {
              lhs[i][j][k][0][m][n] = 0.0;
              lhs[i][j][k][1][m][n] = 0.0;
              lhs[i][j][k][2][m][n] = 0.0;
            }
          }
        }
      }
    }
    
#pragma omp parallel for private (i,j,k,m)
    for (i = 0; i <= grid_points[0] - 1; i += 1) {
      
#pragma omp parallel for private (j,k,m)
      for (j = 0; j <= grid_points[1] - 1; j += 1) {
        
#pragma omp parallel for private (k,m)
        for (k = 0; k <= grid_points[2] - 1; k += 1) {
          
#pragma omp parallel for private (m)
          for (m = 0; m <= 4; m += 1) {
            lhs[i][j][k][1][m][m] = 1.0;
          }
        }
      }
    }
  }
}

static void lhsx()
{
  int i;
  int j;
  int k;
  
#pragma omp parallel for private (tmp1,tmp2,tmp3,i,j,k)
  for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
    
#pragma omp parallel for private (tmp1,tmp2,tmp3,i,k)
    for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
      
#pragma omp parallel for private (tmp1,tmp2,tmp3,i) firstprivate (c3c4,c1345,c1,c2,con43)
      for (i = 0; i <= grid_points[0] - 1; i += 1) {
        tmp1 = 1.0 / u[i][j][k][0];
        tmp2 = tmp1 * tmp1;
        tmp3 = tmp1 * tmp2;
        fjac[i][j][k][0][0] = 0.0;
        fjac[i][j][k][0][1] = 1.0;
        fjac[i][j][k][0][2] = 0.0;
        fjac[i][j][k][0][3] = 0.0;
        fjac[i][j][k][0][4] = 0.0;
        fjac[i][j][k][1][0] = -(u[i][j][k][1] * tmp2 * u[i][j][k][1]) + c2 * 0.50 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2;
        fjac[i][j][k][1][1] = (2.0 - c2) * (u[i][j][k][1] / u[i][j][k][0]);
        fjac[i][j][k][1][2] = -c2 * (u[i][j][k][2] * tmp1);
        fjac[i][j][k][1][3] = -c2 * (u[i][j][k][3] * tmp1);
        fjac[i][j][k][1][4] = c2;
        fjac[i][j][k][2][0] = -(u[i][j][k][1] * u[i][j][k][2]) * tmp2;
        fjac[i][j][k][2][1] = u[i][j][k][2] * tmp1;
        fjac[i][j][k][2][2] = u[i][j][k][1] * tmp1;
        fjac[i][j][k][2][3] = 0.0;
        fjac[i][j][k][2][4] = 0.0;
        fjac[i][j][k][3][0] = -(u[i][j][k][1] * u[i][j][k][3]) * tmp2;
        fjac[i][j][k][3][1] = u[i][j][k][3] * tmp1;
        fjac[i][j][k][3][2] = 0.0;
        fjac[i][j][k][3][3] = u[i][j][k][1] * tmp1;
        fjac[i][j][k][3][4] = 0.0;
        fjac[i][j][k][4][0] = (c2 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2 - c1 * (u[i][j][k][4] * tmp1)) * (u[i][j][k][1] * tmp1);
        fjac[i][j][k][4][1] = c1 * u[i][j][k][4] * tmp1 - 0.50 * c2 * (3.0 * u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2;
        fjac[i][j][k][4][2] = -c2 * (u[i][j][k][2] * u[i][j][k][1]) * tmp2;
        fjac[i][j][k][4][3] = -c2 * (u[i][j][k][3] * u[i][j][k][1]) * tmp2;
        fjac[i][j][k][4][4] = c1 * (u[i][j][k][1] * tmp1);
        njac[i][j][k][0][0] = 0.0;
        njac[i][j][k][0][1] = 0.0;
        njac[i][j][k][0][2] = 0.0;
        njac[i][j][k][0][3] = 0.0;
        njac[i][j][k][0][4] = 0.0;
        njac[i][j][k][1][0] = -con43 * c3c4 * tmp2 * u[i][j][k][1];
        njac[i][j][k][1][1] = con43 * c3c4 * tmp1;
        njac[i][j][k][1][2] = 0.0;
        njac[i][j][k][1][3] = 0.0;
        njac[i][j][k][1][4] = 0.0;
        njac[i][j][k][2][0] = -c3c4 * tmp2 * u[i][j][k][2];
        njac[i][j][k][2][1] = 0.0;
        njac[i][j][k][2][2] = c3c4 * tmp1;
        njac[i][j][k][2][3] = 0.0;
        njac[i][j][k][2][4] = 0.0;
        njac[i][j][k][3][0] = -c3c4 * tmp2 * u[i][j][k][3];
        njac[i][j][k][3][1] = 0.0;
        njac[i][j][k][3][2] = 0.0;
        njac[i][j][k][3][3] = c3c4 * tmp1;
        njac[i][j][k][3][4] = 0.0;
        njac[i][j][k][4][0] = -(con43 * c3c4 - c1345) * tmp3 * (u[i][j][k][1] * u[i][j][k][1]) - (c3c4 - c1345) * tmp3 * (u[i][j][k][2] * u[i][j][k][2]) - (c3c4 - c1345) * tmp3 * (u[i][j][k][3] * u[i][j][k][3]) - c1345 * tmp2 * u[i][j][k][4];
        njac[i][j][k][4][1] = (con43 * c3c4 - c1345) * tmp2 * u[i][j][k][1];
        njac[i][j][k][4][2] = (c3c4 - c1345) * tmp2 * u[i][j][k][2];
        njac[i][j][k][4][3] = (c3c4 - c1345) * tmp2 * u[i][j][k][3];
        njac[i][j][k][4][4] = c1345 * tmp1;
      }
      
#pragma omp parallel for private (tmp1,tmp2,i) firstprivate (tx1,tx2,dx1,dx2,dx3,dx4,dx5,dt)
      for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
        tmp1 = dt * tx1;
        tmp2 = dt * tx2;
        lhs[i][j][k][0][0][0] = -tmp2 * fjac[i - 1][j][k][0][0] - tmp1 * njac[i - 1][j][k][0][0] - tmp1 * dx1;
        lhs[i][j][k][0][0][1] = -tmp2 * fjac[i - 1][j][k][0][1] - tmp1 * njac[i - 1][j][k][0][1];
        lhs[i][j][k][0][0][2] = -tmp2 * fjac[i - 1][j][k][0][2] - tmp1 * njac[i - 1][j][k][0][2];
        lhs[i][j][k][0][0][3] = -tmp2 * fjac[i - 1][j][k][0][3] - tmp1 * njac[i - 1][j][k][0][3];
        lhs[i][j][k][0][0][4] = -tmp2 * fjac[i - 1][j][k][0][4] - tmp1 * njac[i - 1][j][k][0][4];
        lhs[i][j][k][0][1][0] = -tmp2 * fjac[i - 1][j][k][1][0] - tmp1 * njac[i - 1][j][k][1][0];
        lhs[i][j][k][0][1][1] = -tmp2 * fjac[i - 1][j][k][1][1] - tmp1 * njac[i - 1][j][k][1][1] - tmp1 * dx2;
        lhs[i][j][k][0][1][2] = -tmp2 * fjac[i - 1][j][k][1][2] - tmp1 * njac[i - 1][j][k][1][2];
        lhs[i][j][k][0][1][3] = -tmp2 * fjac[i - 1][j][k][1][3] - tmp1 * njac[i - 1][j][k][1][3];
        lhs[i][j][k][0][1][4] = -tmp2 * fjac[i - 1][j][k][1][4] - tmp1 * njac[i - 1][j][k][1][4];
        lhs[i][j][k][0][2][0] = -tmp2 * fjac[i - 1][j][k][2][0] - tmp1 * njac[i - 1][j][k][2][0];
        lhs[i][j][k][0][2][1] = -tmp2 * fjac[i - 1][j][k][2][1] - tmp1 * njac[i - 1][j][k][2][1];
        lhs[i][j][k][0][2][2] = -tmp2 * fjac[i - 1][j][k][2][2] - tmp1 * njac[i - 1][j][k][2][2] - tmp1 * dx3;
        lhs[i][j][k][0][2][3] = -tmp2 * fjac[i - 1][j][k][2][3] - tmp1 * njac[i - 1][j][k][2][3];
        lhs[i][j][k][0][2][4] = -tmp2 * fjac[i - 1][j][k][2][4] - tmp1 * njac[i - 1][j][k][2][4];
        lhs[i][j][k][0][3][0] = -tmp2 * fjac[i - 1][j][k][3][0] - tmp1 * njac[i - 1][j][k][3][0];
        lhs[i][j][k][0][3][1] = -tmp2 * fjac[i - 1][j][k][3][1] - tmp1 * njac[i - 1][j][k][3][1];
        lhs[i][j][k][0][3][2] = -tmp2 * fjac[i - 1][j][k][3][2] - tmp1 * njac[i - 1][j][k][3][2];
        lhs[i][j][k][0][3][3] = -tmp2 * fjac[i - 1][j][k][3][3] - tmp1 * njac[i - 1][j][k][3][3] - tmp1 * dx4;
        lhs[i][j][k][0][3][4] = -tmp2 * fjac[i - 1][j][k][3][4] - tmp1 * njac[i - 1][j][k][3][4];
        lhs[i][j][k][0][4][0] = -tmp2 * fjac[i - 1][j][k][4][0] - tmp1 * njac[i - 1][j][k][4][0];
        lhs[i][j][k][0][4][1] = -tmp2 * fjac[i - 1][j][k][4][1] - tmp1 * njac[i - 1][j][k][4][1];
        lhs[i][j][k][0][4][2] = -tmp2 * fjac[i - 1][j][k][4][2] - tmp1 * njac[i - 1][j][k][4][2];
        lhs[i][j][k][0][4][3] = -tmp2 * fjac[i - 1][j][k][4][3] - tmp1 * njac[i - 1][j][k][4][3];
        lhs[i][j][k][0][4][4] = -tmp2 * fjac[i - 1][j][k][4][4] - tmp1 * njac[i - 1][j][k][4][4] - tmp1 * dx5;
        lhs[i][j][k][1][0][0] = 1.0 + tmp1 * 2.0 * njac[i][j][k][0][0] + tmp1 * 2.0 * dx1;
        lhs[i][j][k][1][0][1] = tmp1 * 2.0 * njac[i][j][k][0][1];
        lhs[i][j][k][1][0][2] = tmp1 * 2.0 * njac[i][j][k][0][2];
        lhs[i][j][k][1][0][3] = tmp1 * 2.0 * njac[i][j][k][0][3];
        lhs[i][j][k][1][0][4] = tmp1 * 2.0 * njac[i][j][k][0][4];
        lhs[i][j][k][1][1][0] = tmp1 * 2.0 * njac[i][j][k][1][0];
        lhs[i][j][k][1][1][1] = 1.0 + tmp1 * 2.0 * njac[i][j][k][1][1] + tmp1 * 2.0 * dx2;
        lhs[i][j][k][1][1][2] = tmp1 * 2.0 * njac[i][j][k][1][2];
        lhs[i][j][k][1][1][3] = tmp1 * 2.0 * njac[i][j][k][1][3];
        lhs[i][j][k][1][1][4] = tmp1 * 2.0 * njac[i][j][k][1][4];
        lhs[i][j][k][1][2][0] = tmp1 * 2.0 * njac[i][j][k][2][0];
        lhs[i][j][k][1][2][1] = tmp1 * 2.0 * njac[i][j][k][2][1];
        lhs[i][j][k][1][2][2] = 1.0 + tmp1 * 2.0 * njac[i][j][k][2][2] + tmp1 * 2.0 * dx3;
        lhs[i][j][k][1][2][3] = tmp1 * 2.0 * njac[i][j][k][2][3];
        lhs[i][j][k][1][2][4] = tmp1 * 2.0 * njac[i][j][k][2][4];
        lhs[i][j][k][1][3][0] = tmp1 * 2.0 * njac[i][j][k][3][0];
        lhs[i][j][k][1][3][1] = tmp1 * 2.0 * njac[i][j][k][3][1];
        lhs[i][j][k][1][3][2] = tmp1 * 2.0 * njac[i][j][k][3][2];
        lhs[i][j][k][1][3][3] = 1.0 + tmp1 * 2.0 * njac[i][j][k][3][3] + tmp1 * 2.0 * dx4;
        lhs[i][j][k][1][3][4] = tmp1 * 2.0 * njac[i][j][k][3][4];
        lhs[i][j][k][1][4][0] = tmp1 * 2.0 * njac[i][j][k][4][0];
        lhs[i][j][k][1][4][1] = tmp1 * 2.0 * njac[i][j][k][4][1];
        lhs[i][j][k][1][4][2] = tmp1 * 2.0 * njac[i][j][k][4][2];
        lhs[i][j][k][1][4][3] = tmp1 * 2.0 * njac[i][j][k][4][3];
        lhs[i][j][k][1][4][4] = 1.0 + tmp1 * 2.0 * njac[i][j][k][4][4] + tmp1 * 2.0 * dx5;
        lhs[i][j][k][2][0][0] = tmp2 * fjac[i + 1][j][k][0][0] - tmp1 * njac[i + 1][j][k][0][0] - tmp1 * dx1;
        lhs[i][j][k][2][0][1] = tmp2 * fjac[i + 1][j][k][0][1] - tmp1 * njac[i + 1][j][k][0][1];
        lhs[i][j][k][2][0][2] = tmp2 * fjac[i + 1][j][k][0][2] - tmp1 * njac[i + 1][j][k][0][2];
        lhs[i][j][k][2][0][3] = tmp2 * fjac[i + 1][j][k][0][3] - tmp1 * njac[i + 1][j][k][0][3];
        lhs[i][j][k][2][0][4] = tmp2 * fjac[i + 1][j][k][0][4] - tmp1 * njac[i + 1][j][k][0][4];
        lhs[i][j][k][2][1][0] = tmp2 * fjac[i + 1][j][k][1][0] - tmp1 * njac[i + 1][j][k][1][0];
        lhs[i][j][k][2][1][1] = tmp2 * fjac[i + 1][j][k][1][1] - tmp1 * njac[i + 1][j][k][1][1] - tmp1 * dx2;
        lhs[i][j][k][2][1][2] = tmp2 * fjac[i + 1][j][k][1][2] - tmp1 * njac[i + 1][j][k][1][2];
        lhs[i][j][k][2][1][3] = tmp2 * fjac[i + 1][j][k][1][3] - tmp1 * njac[i + 1][j][k][1][3];
        lhs[i][j][k][2][1][4] = tmp2 * fjac[i + 1][j][k][1][4] - tmp1 * njac[i + 1][j][k][1][4];
        lhs[i][j][k][2][2][0] = tmp2 * fjac[i + 1][j][k][2][0] - tmp1 * njac[i + 1][j][k][2][0];
        lhs[i][j][k][2][2][1] = tmp2 * fjac[i + 1][j][k][2][1] - tmp1 * njac[i + 1][j][k][2][1];
        lhs[i][j][k][2][2][2] = tmp2 * fjac[i + 1][j][k][2][2] - tmp1 * njac[i + 1][j][k][2][2] - tmp1 * dx3;
        lhs[i][j][k][2][2][3] = tmp2 * fjac[i + 1][j][k][2][3] - tmp1 * njac[i + 1][j][k][2][3];
        lhs[i][j][k][2][2][4] = tmp2 * fjac[i + 1][j][k][2][4] - tmp1 * njac[i + 1][j][k][2][4];
        lhs[i][j][k][2][3][0] = tmp2 * fjac[i + 1][j][k][3][0] - tmp1 * njac[i + 1][j][k][3][0];
        lhs[i][j][k][2][3][1] = tmp2 * fjac[i + 1][j][k][3][1] - tmp1 * njac[i + 1][j][k][3][1];
        lhs[i][j][k][2][3][2] = tmp2 * fjac[i + 1][j][k][3][2] - tmp1 * njac[i + 1][j][k][3][2];
        lhs[i][j][k][2][3][3] = tmp2 * fjac[i + 1][j][k][3][3] - tmp1 * njac[i + 1][j][k][3][3] - tmp1 * dx4;
        lhs[i][j][k][2][3][4] = tmp2 * fjac[i + 1][j][k][3][4] - tmp1 * njac[i + 1][j][k][3][4];
        lhs[i][j][k][2][4][0] = tmp2 * fjac[i + 1][j][k][4][0] - tmp1 * njac[i + 1][j][k][4][0];
        lhs[i][j][k][2][4][1] = tmp2 * fjac[i + 1][j][k][4][1] - tmp1 * njac[i + 1][j][k][4][1];
        lhs[i][j][k][2][4][2] = tmp2 * fjac[i + 1][j][k][4][2] - tmp1 * njac[i + 1][j][k][4][2];
        lhs[i][j][k][2][4][3] = tmp2 * fjac[i + 1][j][k][4][3] - tmp1 * njac[i + 1][j][k][4][3];
        lhs[i][j][k][2][4][4] = tmp2 * fjac[i + 1][j][k][4][4] - tmp1 * njac[i + 1][j][k][4][4] - tmp1 * dx5;
      }
    }
  }
}

static void lhsy()
{
  int i;
  int j;
  int k;
  
#pragma omp parallel for private (tmp1,tmp2,tmp3,i,j,k)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (tmp1,tmp2,tmp3,j,k)
    for (j = 0; j <= grid_points[1] - 1; j += 1) {
      
#pragma omp parallel for private (tmp1,tmp2,tmp3,k) firstprivate (c3c4,c1345,c1,c2,con43)
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        tmp1 = 1.0 / u[i][j][k][0];
        tmp2 = tmp1 * tmp1;
        tmp3 = tmp1 * tmp2;
        fjac[i][j][k][0][0] = 0.0;
        fjac[i][j][k][0][1] = 0.0;
        fjac[i][j][k][0][2] = 1.0;
        fjac[i][j][k][0][3] = 0.0;
        fjac[i][j][k][0][4] = 0.0;
        fjac[i][j][k][1][0] = -(u[i][j][k][1] * u[i][j][k][2]) * tmp2;
        fjac[i][j][k][1][1] = u[i][j][k][2] * tmp1;
        fjac[i][j][k][1][2] = u[i][j][k][1] * tmp1;
        fjac[i][j][k][1][3] = 0.0;
        fjac[i][j][k][1][4] = 0.0;
        fjac[i][j][k][2][0] = -(u[i][j][k][2] * u[i][j][k][2] * tmp2) + 0.50 * c2 * ((u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2);
        fjac[i][j][k][2][1] = -c2 * u[i][j][k][1] * tmp1;
        fjac[i][j][k][2][2] = (2.0 - c2) * u[i][j][k][2] * tmp1;
        fjac[i][j][k][2][3] = -c2 * u[i][j][k][3] * tmp1;
        fjac[i][j][k][2][4] = c2;
        fjac[i][j][k][3][0] = -(u[i][j][k][2] * u[i][j][k][3]) * tmp2;
        fjac[i][j][k][3][1] = 0.0;
        fjac[i][j][k][3][2] = u[i][j][k][3] * tmp1;
        fjac[i][j][k][3][3] = u[i][j][k][2] * tmp1;
        fjac[i][j][k][3][4] = 0.0;
        fjac[i][j][k][4][0] = (c2 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2 - c1 * u[i][j][k][4] * tmp1) * u[i][j][k][2] * tmp1;
        fjac[i][j][k][4][1] = -c2 * u[i][j][k][1] * u[i][j][k][2] * tmp2;
        fjac[i][j][k][4][2] = c1 * u[i][j][k][4] * tmp1 - 0.50 * c2 * ((u[i][j][k][1] * u[i][j][k][1] + 3.0 * u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2);
        fjac[i][j][k][4][3] = -c2 * (u[i][j][k][2] * u[i][j][k][3]) * tmp2;
        fjac[i][j][k][4][4] = c1 * u[i][j][k][2] * tmp1;
        njac[i][j][k][0][0] = 0.0;
        njac[i][j][k][0][1] = 0.0;
        njac[i][j][k][0][2] = 0.0;
        njac[i][j][k][0][3] = 0.0;
        njac[i][j][k][0][4] = 0.0;
        njac[i][j][k][1][0] = -c3c4 * tmp2 * u[i][j][k][1];
        njac[i][j][k][1][1] = c3c4 * tmp1;
        njac[i][j][k][1][2] = 0.0;
        njac[i][j][k][1][3] = 0.0;
        njac[i][j][k][1][4] = 0.0;
        njac[i][j][k][2][0] = -con43 * c3c4 * tmp2 * u[i][j][k][2];
        njac[i][j][k][2][1] = 0.0;
        njac[i][j][k][2][2] = con43 * c3c4 * tmp1;
        njac[i][j][k][2][3] = 0.0;
        njac[i][j][k][2][4] = 0.0;
        njac[i][j][k][3][0] = -c3c4 * tmp2 * u[i][j][k][3];
        njac[i][j][k][3][1] = 0.0;
        njac[i][j][k][3][2] = 0.0;
        njac[i][j][k][3][3] = c3c4 * tmp1;
        njac[i][j][k][3][4] = 0.0;
        njac[i][j][k][4][0] = -(c3c4 - c1345) * tmp3 * (u[i][j][k][1] * u[i][j][k][1]) - (con43 * c3c4 - c1345) * tmp3 * (u[i][j][k][2] * u[i][j][k][2]) - (c3c4 - c1345) * tmp3 * (u[i][j][k][3] * u[i][j][k][3]) - c1345 * tmp2 * u[i][j][k][4];
        njac[i][j][k][4][1] = (c3c4 - c1345) * tmp2 * u[i][j][k][1];
        njac[i][j][k][4][2] = (con43 * c3c4 - c1345) * tmp2 * u[i][j][k][2];
        njac[i][j][k][4][3] = (c3c4 - c1345) * tmp2 * u[i][j][k][3];
        njac[i][j][k][4][4] = c1345 * tmp1;
      }
    }
  }
  
#pragma omp parallel for private (tmp1,tmp2,i,j,k)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (tmp1,tmp2,j,k)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      
#pragma omp parallel for private (tmp1,tmp2,k) firstprivate (ty1,ty2,dy1,dy2,dy3,dy4,dy5,dt)
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        tmp1 = dt * ty1;
        tmp2 = dt * ty2;
        lhs[i][j][k][0][0][0] = -tmp2 * fjac[i][j - 1][k][0][0] - tmp1 * njac[i][j - 1][k][0][0] - tmp1 * dy1;
        lhs[i][j][k][0][0][1] = -tmp2 * fjac[i][j - 1][k][0][1] - tmp1 * njac[i][j - 1][k][0][1];
        lhs[i][j][k][0][0][2] = -tmp2 * fjac[i][j - 1][k][0][2] - tmp1 * njac[i][j - 1][k][0][2];
        lhs[i][j][k][0][0][3] = -tmp2 * fjac[i][j - 1][k][0][3] - tmp1 * njac[i][j - 1][k][0][3];
        lhs[i][j][k][0][0][4] = -tmp2 * fjac[i][j - 1][k][0][4] - tmp1 * njac[i][j - 1][k][0][4];
        lhs[i][j][k][0][1][0] = -tmp2 * fjac[i][j - 1][k][1][0] - tmp1 * njac[i][j - 1][k][1][0];
        lhs[i][j][k][0][1][1] = -tmp2 * fjac[i][j - 1][k][1][1] - tmp1 * njac[i][j - 1][k][1][1] - tmp1 * dy2;
        lhs[i][j][k][0][1][2] = -tmp2 * fjac[i][j - 1][k][1][2] - tmp1 * njac[i][j - 1][k][1][2];
        lhs[i][j][k][0][1][3] = -tmp2 * fjac[i][j - 1][k][1][3] - tmp1 * njac[i][j - 1][k][1][3];
        lhs[i][j][k][0][1][4] = -tmp2 * fjac[i][j - 1][k][1][4] - tmp1 * njac[i][j - 1][k][1][4];
        lhs[i][j][k][0][2][0] = -tmp2 * fjac[i][j - 1][k][2][0] - tmp1 * njac[i][j - 1][k][2][0];
        lhs[i][j][k][0][2][1] = -tmp2 * fjac[i][j - 1][k][2][1] - tmp1 * njac[i][j - 1][k][2][1];
        lhs[i][j][k][0][2][2] = -tmp2 * fjac[i][j - 1][k][2][2] - tmp1 * njac[i][j - 1][k][2][2] - tmp1 * dy3;
        lhs[i][j][k][0][2][3] = -tmp2 * fjac[i][j - 1][k][2][3] - tmp1 * njac[i][j - 1][k][2][3];
        lhs[i][j][k][0][2][4] = -tmp2 * fjac[i][j - 1][k][2][4] - tmp1 * njac[i][j - 1][k][2][4];
        lhs[i][j][k][0][3][0] = -tmp2 * fjac[i][j - 1][k][3][0] - tmp1 * njac[i][j - 1][k][3][0];
        lhs[i][j][k][0][3][1] = -tmp2 * fjac[i][j - 1][k][3][1] - tmp1 * njac[i][j - 1][k][3][1];
        lhs[i][j][k][0][3][2] = -tmp2 * fjac[i][j - 1][k][3][2] - tmp1 * njac[i][j - 1][k][3][2];
        lhs[i][j][k][0][3][3] = -tmp2 * fjac[i][j - 1][k][3][3] - tmp1 * njac[i][j - 1][k][3][3] - tmp1 * dy4;
        lhs[i][j][k][0][3][4] = -tmp2 * fjac[i][j - 1][k][3][4] - tmp1 * njac[i][j - 1][k][3][4];
        lhs[i][j][k][0][4][0] = -tmp2 * fjac[i][j - 1][k][4][0] - tmp1 * njac[i][j - 1][k][4][0];
        lhs[i][j][k][0][4][1] = -tmp2 * fjac[i][j - 1][k][4][1] - tmp1 * njac[i][j - 1][k][4][1];
        lhs[i][j][k][0][4][2] = -tmp2 * fjac[i][j - 1][k][4][2] - tmp1 * njac[i][j - 1][k][4][2];
        lhs[i][j][k][0][4][3] = -tmp2 * fjac[i][j - 1][k][4][3] - tmp1 * njac[i][j - 1][k][4][3];
        lhs[i][j][k][0][4][4] = -tmp2 * fjac[i][j - 1][k][4][4] - tmp1 * njac[i][j - 1][k][4][4] - tmp1 * dy5;
        lhs[i][j][k][1][0][0] = 1.0 + tmp1 * 2.0 * njac[i][j][k][0][0] + tmp1 * 2.0 * dy1;
        lhs[i][j][k][1][0][1] = tmp1 * 2.0 * njac[i][j][k][0][1];
        lhs[i][j][k][1][0][2] = tmp1 * 2.0 * njac[i][j][k][0][2];
        lhs[i][j][k][1][0][3] = tmp1 * 2.0 * njac[i][j][k][0][3];
        lhs[i][j][k][1][0][4] = tmp1 * 2.0 * njac[i][j][k][0][4];
        lhs[i][j][k][1][1][0] = tmp1 * 2.0 * njac[i][j][k][1][0];
        lhs[i][j][k][1][1][1] = 1.0 + tmp1 * 2.0 * njac[i][j][k][1][1] + tmp1 * 2.0 * dy2;
        lhs[i][j][k][1][1][2] = tmp1 * 2.0 * njac[i][j][k][1][2];
        lhs[i][j][k][1][1][3] = tmp1 * 2.0 * njac[i][j][k][1][3];
        lhs[i][j][k][1][1][4] = tmp1 * 2.0 * njac[i][j][k][1][4];
        lhs[i][j][k][1][2][0] = tmp1 * 2.0 * njac[i][j][k][2][0];
        lhs[i][j][k][1][2][1] = tmp1 * 2.0 * njac[i][j][k][2][1];
        lhs[i][j][k][1][2][2] = 1.0 + tmp1 * 2.0 * njac[i][j][k][2][2] + tmp1 * 2.0 * dy3;
        lhs[i][j][k][1][2][3] = tmp1 * 2.0 * njac[i][j][k][2][3];
        lhs[i][j][k][1][2][4] = tmp1 * 2.0 * njac[i][j][k][2][4];
        lhs[i][j][k][1][3][0] = tmp1 * 2.0 * njac[i][j][k][3][0];
        lhs[i][j][k][1][3][1] = tmp1 * 2.0 * njac[i][j][k][3][1];
        lhs[i][j][k][1][3][2] = tmp1 * 2.0 * njac[i][j][k][3][2];
        lhs[i][j][k][1][3][3] = 1.0 + tmp1 * 2.0 * njac[i][j][k][3][3] + tmp1 * 2.0 * dy4;
        lhs[i][j][k][1][3][4] = tmp1 * 2.0 * njac[i][j][k][3][4];
        lhs[i][j][k][1][4][0] = tmp1 * 2.0 * njac[i][j][k][4][0];
        lhs[i][j][k][1][4][1] = tmp1 * 2.0 * njac[i][j][k][4][1];
        lhs[i][j][k][1][4][2] = tmp1 * 2.0 * njac[i][j][k][4][2];
        lhs[i][j][k][1][4][3] = tmp1 * 2.0 * njac[i][j][k][4][3];
        lhs[i][j][k][1][4][4] = 1.0 + tmp1 * 2.0 * njac[i][j][k][4][4] + tmp1 * 2.0 * dy5;
        lhs[i][j][k][2][0][0] = tmp2 * fjac[i][j + 1][k][0][0] - tmp1 * njac[i][j + 1][k][0][0] - tmp1 * dy1;
        lhs[i][j][k][2][0][1] = tmp2 * fjac[i][j + 1][k][0][1] - tmp1 * njac[i][j + 1][k][0][1];
        lhs[i][j][k][2][0][2] = tmp2 * fjac[i][j + 1][k][0][2] - tmp1 * njac[i][j + 1][k][0][2];
        lhs[i][j][k][2][0][3] = tmp2 * fjac[i][j + 1][k][0][3] - tmp1 * njac[i][j + 1][k][0][3];
        lhs[i][j][k][2][0][4] = tmp2 * fjac[i][j + 1][k][0][4] - tmp1 * njac[i][j + 1][k][0][4];
        lhs[i][j][k][2][1][0] = tmp2 * fjac[i][j + 1][k][1][0] - tmp1 * njac[i][j + 1][k][1][0];
        lhs[i][j][k][2][1][1] = tmp2 * fjac[i][j + 1][k][1][1] - tmp1 * njac[i][j + 1][k][1][1] - tmp1 * dy2;
        lhs[i][j][k][2][1][2] = tmp2 * fjac[i][j + 1][k][1][2] - tmp1 * njac[i][j + 1][k][1][2];
        lhs[i][j][k][2][1][3] = tmp2 * fjac[i][j + 1][k][1][3] - tmp1 * njac[i][j + 1][k][1][3];
        lhs[i][j][k][2][1][4] = tmp2 * fjac[i][j + 1][k][1][4] - tmp1 * njac[i][j + 1][k][1][4];
        lhs[i][j][k][2][2][0] = tmp2 * fjac[i][j + 1][k][2][0] - tmp1 * njac[i][j + 1][k][2][0];
        lhs[i][j][k][2][2][1] = tmp2 * fjac[i][j + 1][k][2][1] - tmp1 * njac[i][j + 1][k][2][1];
        lhs[i][j][k][2][2][2] = tmp2 * fjac[i][j + 1][k][2][2] - tmp1 * njac[i][j + 1][k][2][2] - tmp1 * dy3;
        lhs[i][j][k][2][2][3] = tmp2 * fjac[i][j + 1][k][2][3] - tmp1 * njac[i][j + 1][k][2][3];
        lhs[i][j][k][2][2][4] = tmp2 * fjac[i][j + 1][k][2][4] - tmp1 * njac[i][j + 1][k][2][4];
        lhs[i][j][k][2][3][0] = tmp2 * fjac[i][j + 1][k][3][0] - tmp1 * njac[i][j + 1][k][3][0];
        lhs[i][j][k][2][3][1] = tmp2 * fjac[i][j + 1][k][3][1] - tmp1 * njac[i][j + 1][k][3][1];
        lhs[i][j][k][2][3][2] = tmp2 * fjac[i][j + 1][k][3][2] - tmp1 * njac[i][j + 1][k][3][2];
        lhs[i][j][k][2][3][3] = tmp2 * fjac[i][j + 1][k][3][3] - tmp1 * njac[i][j + 1][k][3][3] - tmp1 * dy4;
        lhs[i][j][k][2][3][4] = tmp2 * fjac[i][j + 1][k][3][4] - tmp1 * njac[i][j + 1][k][3][4];
        lhs[i][j][k][2][4][0] = tmp2 * fjac[i][j + 1][k][4][0] - tmp1 * njac[i][j + 1][k][4][0];
        lhs[i][j][k][2][4][1] = tmp2 * fjac[i][j + 1][k][4][1] - tmp1 * njac[i][j + 1][k][4][1];
        lhs[i][j][k][2][4][2] = tmp2 * fjac[i][j + 1][k][4][2] - tmp1 * njac[i][j + 1][k][4][2];
        lhs[i][j][k][2][4][3] = tmp2 * fjac[i][j + 1][k][4][3] - tmp1 * njac[i][j + 1][k][4][3];
        lhs[i][j][k][2][4][4] = tmp2 * fjac[i][j + 1][k][4][4] - tmp1 * njac[i][j + 1][k][4][4] - tmp1 * dy5;
      }
    }
  }
}

static void lhsz()
{
  int i;
  int j;
  int k;
  
#pragma omp parallel for private (tmp1,tmp2,tmp3,i,j,k)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (tmp1,tmp2,tmp3,j,k)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      
#pragma omp parallel for private (tmp1,tmp2,tmp3,k) firstprivate (c3c4,c1345,c1,c2,c3,c4,con43)
      for (k = 0; k <= grid_points[2] - 1; k += 1) {
        tmp1 = 1.0 / u[i][j][k][0];
        tmp2 = tmp1 * tmp1;
        tmp3 = tmp1 * tmp2;
        fjac[i][j][k][0][0] = 0.0;
        fjac[i][j][k][0][1] = 0.0;
        fjac[i][j][k][0][2] = 0.0;
        fjac[i][j][k][0][3] = 1.0;
        fjac[i][j][k][0][4] = 0.0;
        fjac[i][j][k][1][0] = -(u[i][j][k][1] * u[i][j][k][3]) * tmp2;
        fjac[i][j][k][1][1] = u[i][j][k][3] * tmp1;
        fjac[i][j][k][1][2] = 0.0;
        fjac[i][j][k][1][3] = u[i][j][k][1] * tmp1;
        fjac[i][j][k][1][4] = 0.0;
        fjac[i][j][k][2][0] = -(u[i][j][k][2] * u[i][j][k][3]) * tmp2;
        fjac[i][j][k][2][1] = 0.0;
        fjac[i][j][k][2][2] = u[i][j][k][3] * tmp1;
        fjac[i][j][k][2][3] = u[i][j][k][2] * tmp1;
        fjac[i][j][k][2][4] = 0.0;
        fjac[i][j][k][3][0] = -(u[i][j][k][3] * u[i][j][k][3] * tmp2) + 0.50 * c2 * ((u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2);
        fjac[i][j][k][3][1] = -c2 * u[i][j][k][1] * tmp1;
        fjac[i][j][k][3][2] = -c2 * u[i][j][k][2] * tmp1;
        fjac[i][j][k][3][3] = (2.0 - c2) * u[i][j][k][3] * tmp1;
        fjac[i][j][k][3][4] = c2;
        fjac[i][j][k][4][0] = (c2 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2 - c1 * (u[i][j][k][4] * tmp1)) * (u[i][j][k][3] * tmp1);
        fjac[i][j][k][4][1] = -c2 * (u[i][j][k][1] * u[i][j][k][3]) * tmp2;
        fjac[i][j][k][4][2] = -c2 * (u[i][j][k][2] * u[i][j][k][3]) * tmp2;
        fjac[i][j][k][4][3] = c1 * (u[i][j][k][4] * tmp1) - 0.50 * c2 * ((u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + 3.0 * u[i][j][k][3] * u[i][j][k][3]) * tmp2);
        fjac[i][j][k][4][4] = c1 * u[i][j][k][3] * tmp1;
        njac[i][j][k][0][0] = 0.0;
        njac[i][j][k][0][1] = 0.0;
        njac[i][j][k][0][2] = 0.0;
        njac[i][j][k][0][3] = 0.0;
        njac[i][j][k][0][4] = 0.0;
        njac[i][j][k][1][0] = -c3c4 * tmp2 * u[i][j][k][1];
        njac[i][j][k][1][1] = c3c4 * tmp1;
        njac[i][j][k][1][2] = 0.0;
        njac[i][j][k][1][3] = 0.0;
        njac[i][j][k][1][4] = 0.0;
        njac[i][j][k][2][0] = -c3c4 * tmp2 * u[i][j][k][2];
        njac[i][j][k][2][1] = 0.0;
        njac[i][j][k][2][2] = c3c4 * tmp1;
        njac[i][j][k][2][3] = 0.0;
        njac[i][j][k][2][4] = 0.0;
        njac[i][j][k][3][0] = -con43 * c3c4 * tmp2 * u[i][j][k][3];
        njac[i][j][k][3][1] = 0.0;
        njac[i][j][k][3][2] = 0.0;
        njac[i][j][k][3][3] = con43 * c3 * c4 * tmp1;
        njac[i][j][k][3][4] = 0.0;
        njac[i][j][k][4][0] = -(c3c4 - c1345) * tmp3 * (u[i][j][k][1] * u[i][j][k][1]) - (c3c4 - c1345) * tmp3 * (u[i][j][k][2] * u[i][j][k][2]) - (con43 * c3c4 - c1345) * tmp3 * (u[i][j][k][3] * u[i][j][k][3]) - c1345 * tmp2 * u[i][j][k][4];
        njac[i][j][k][4][1] = (c3c4 - c1345) * tmp2 * u[i][j][k][1];
        njac[i][j][k][4][2] = (c3c4 - c1345) * tmp2 * u[i][j][k][2];
        njac[i][j][k][4][3] = (con43 * c3c4 - c1345) * tmp2 * u[i][j][k][3];
        njac[i][j][k][4][4] = c1345 * tmp1;
      }
    }
  }
  
#pragma omp parallel for private (tmp1,tmp2,i,j,k)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (tmp1,tmp2,j,k)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      
#pragma omp parallel for private (tmp1,tmp2,k) firstprivate (tz1,tz2,dz1,dz2,dz3,dz4,dz5,dt)
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        tmp1 = dt * tz1;
        tmp2 = dt * tz2;
        lhs[i][j][k][0][0][0] = -tmp2 * fjac[i][j][k - 1][0][0] - tmp1 * njac[i][j][k - 1][0][0] - tmp1 * dz1;
        lhs[i][j][k][0][0][1] = -tmp2 * fjac[i][j][k - 1][0][1] - tmp1 * njac[i][j][k - 1][0][1];
        lhs[i][j][k][0][0][2] = -tmp2 * fjac[i][j][k - 1][0][2] - tmp1 * njac[i][j][k - 1][0][2];
        lhs[i][j][k][0][0][3] = -tmp2 * fjac[i][j][k - 1][0][3] - tmp1 * njac[i][j][k - 1][0][3];
        lhs[i][j][k][0][0][4] = -tmp2 * fjac[i][j][k - 1][0][4] - tmp1 * njac[i][j][k - 1][0][4];
        lhs[i][j][k][0][1][0] = -tmp2 * fjac[i][j][k - 1][1][0] - tmp1 * njac[i][j][k - 1][1][0];
        lhs[i][j][k][0][1][1] = -tmp2 * fjac[i][j][k - 1][1][1] - tmp1 * njac[i][j][k - 1][1][1] - tmp1 * dz2;
        lhs[i][j][k][0][1][2] = -tmp2 * fjac[i][j][k - 1][1][2] - tmp1 * njac[i][j][k - 1][1][2];
        lhs[i][j][k][0][1][3] = -tmp2 * fjac[i][j][k - 1][1][3] - tmp1 * njac[i][j][k - 1][1][3];
        lhs[i][j][k][0][1][4] = -tmp2 * fjac[i][j][k - 1][1][4] - tmp1 * njac[i][j][k - 1][1][4];
        lhs[i][j][k][0][2][0] = -tmp2 * fjac[i][j][k - 1][2][0] - tmp1 * njac[i][j][k - 1][2][0];
        lhs[i][j][k][0][2][1] = -tmp2 * fjac[i][j][k - 1][2][1] - tmp1 * njac[i][j][k - 1][2][1];
        lhs[i][j][k][0][2][2] = -tmp2 * fjac[i][j][k - 1][2][2] - tmp1 * njac[i][j][k - 1][2][2] - tmp1 * dz3;
        lhs[i][j][k][0][2][3] = -tmp2 * fjac[i][j][k - 1][2][3] - tmp1 * njac[i][j][k - 1][2][3];
        lhs[i][j][k][0][2][4] = -tmp2 * fjac[i][j][k - 1][2][4] - tmp1 * njac[i][j][k - 1][2][4];
        lhs[i][j][k][0][3][0] = -tmp2 * fjac[i][j][k - 1][3][0] - tmp1 * njac[i][j][k - 1][3][0];
        lhs[i][j][k][0][3][1] = -tmp2 * fjac[i][j][k - 1][3][1] - tmp1 * njac[i][j][k - 1][3][1];
        lhs[i][j][k][0][3][2] = -tmp2 * fjac[i][j][k - 1][3][2] - tmp1 * njac[i][j][k - 1][3][2];
        lhs[i][j][k][0][3][3] = -tmp2 * fjac[i][j][k - 1][3][3] - tmp1 * njac[i][j][k - 1][3][3] - tmp1 * dz4;
        lhs[i][j][k][0][3][4] = -tmp2 * fjac[i][j][k - 1][3][4] - tmp1 * njac[i][j][k - 1][3][4];
        lhs[i][j][k][0][4][0] = -tmp2 * fjac[i][j][k - 1][4][0] - tmp1 * njac[i][j][k - 1][4][0];
        lhs[i][j][k][0][4][1] = -tmp2 * fjac[i][j][k - 1][4][1] - tmp1 * njac[i][j][k - 1][4][1];
        lhs[i][j][k][0][4][2] = -tmp2 * fjac[i][j][k - 1][4][2] - tmp1 * njac[i][j][k - 1][4][2];
        lhs[i][j][k][0][4][3] = -tmp2 * fjac[i][j][k - 1][4][3] - tmp1 * njac[i][j][k - 1][4][3];
        lhs[i][j][k][0][4][4] = -tmp2 * fjac[i][j][k - 1][4][4] - tmp1 * njac[i][j][k - 1][4][4] - tmp1 * dz5;
        lhs[i][j][k][1][0][0] = 1.0 + tmp1 * 2.0 * njac[i][j][k][0][0] + tmp1 * 2.0 * dz1;
        lhs[i][j][k][1][0][1] = tmp1 * 2.0 * njac[i][j][k][0][1];
        lhs[i][j][k][1][0][2] = tmp1 * 2.0 * njac[i][j][k][0][2];
        lhs[i][j][k][1][0][3] = tmp1 * 2.0 * njac[i][j][k][0][3];
        lhs[i][j][k][1][0][4] = tmp1 * 2.0 * njac[i][j][k][0][4];
        lhs[i][j][k][1][1][0] = tmp1 * 2.0 * njac[i][j][k][1][0];
        lhs[i][j][k][1][1][1] = 1.0 + tmp1 * 2.0 * njac[i][j][k][1][1] + tmp1 * 2.0 * dz2;
        lhs[i][j][k][1][1][2] = tmp1 * 2.0 * njac[i][j][k][1][2];
        lhs[i][j][k][1][1][3] = tmp1 * 2.0 * njac[i][j][k][1][3];
        lhs[i][j][k][1][1][4] = tmp1 * 2.0 * njac[i][j][k][1][4];
        lhs[i][j][k][1][2][0] = tmp1 * 2.0 * njac[i][j][k][2][0];
        lhs[i][j][k][1][2][1] = tmp1 * 2.0 * njac[i][j][k][2][1];
        lhs[i][j][k][1][2][2] = 1.0 + tmp1 * 2.0 * njac[i][j][k][2][2] + tmp1 * 2.0 * dz3;
        lhs[i][j][k][1][2][3] = tmp1 * 2.0 * njac[i][j][k][2][3];
        lhs[i][j][k][1][2][4] = tmp1 * 2.0 * njac[i][j][k][2][4];
        lhs[i][j][k][1][3][0] = tmp1 * 2.0 * njac[i][j][k][3][0];
        lhs[i][j][k][1][3][1] = tmp1 * 2.0 * njac[i][j][k][3][1];
        lhs[i][j][k][1][3][2] = tmp1 * 2.0 * njac[i][j][k][3][2];
        lhs[i][j][k][1][3][3] = 1.0 + tmp1 * 2.0 * njac[i][j][k][3][3] + tmp1 * 2.0 * dz4;
        lhs[i][j][k][1][3][4] = tmp1 * 2.0 * njac[i][j][k][3][4];
        lhs[i][j][k][1][4][0] = tmp1 * 2.0 * njac[i][j][k][4][0];
        lhs[i][j][k][1][4][1] = tmp1 * 2.0 * njac[i][j][k][4][1];
        lhs[i][j][k][1][4][2] = tmp1 * 2.0 * njac[i][j][k][4][2];
        lhs[i][j][k][1][4][3] = tmp1 * 2.0 * njac[i][j][k][4][3];
        lhs[i][j][k][1][4][4] = 1.0 + tmp1 * 2.0 * njac[i][j][k][4][4] + tmp1 * 2.0 * dz5;
        lhs[i][j][k][2][0][0] = tmp2 * fjac[i][j][k + 1][0][0] - tmp1 * njac[i][j][k + 1][0][0] - tmp1 * dz1;
        lhs[i][j][k][2][0][1] = tmp2 * fjac[i][j][k + 1][0][1] - tmp1 * njac[i][j][k + 1][0][1];
        lhs[i][j][k][2][0][2] = tmp2 * fjac[i][j][k + 1][0][2] - tmp1 * njac[i][j][k + 1][0][2];
        lhs[i][j][k][2][0][3] = tmp2 * fjac[i][j][k + 1][0][3] - tmp1 * njac[i][j][k + 1][0][3];
        lhs[i][j][k][2][0][4] = tmp2 * fjac[i][j][k + 1][0][4] - tmp1 * njac[i][j][k + 1][0][4];
        lhs[i][j][k][2][1][0] = tmp2 * fjac[i][j][k + 1][1][0] - tmp1 * njac[i][j][k + 1][1][0];
        lhs[i][j][k][2][1][1] = tmp2 * fjac[i][j][k + 1][1][1] - tmp1 * njac[i][j][k + 1][1][1] - tmp1 * dz2;
        lhs[i][j][k][2][1][2] = tmp2 * fjac[i][j][k + 1][1][2] - tmp1 * njac[i][j][k + 1][1][2];
        lhs[i][j][k][2][1][3] = tmp2 * fjac[i][j][k + 1][1][3] - tmp1 * njac[i][j][k + 1][1][3];
        lhs[i][j][k][2][1][4] = tmp2 * fjac[i][j][k + 1][1][4] - tmp1 * njac[i][j][k + 1][1][4];
        lhs[i][j][k][2][2][0] = tmp2 * fjac[i][j][k + 1][2][0] - tmp1 * njac[i][j][k + 1][2][0];
        lhs[i][j][k][2][2][1] = tmp2 * fjac[i][j][k + 1][2][1] - tmp1 * njac[i][j][k + 1][2][1];
        lhs[i][j][k][2][2][2] = tmp2 * fjac[i][j][k + 1][2][2] - tmp1 * njac[i][j][k + 1][2][2] - tmp1 * dz3;
        lhs[i][j][k][2][2][3] = tmp2 * fjac[i][j][k + 1][2][3] - tmp1 * njac[i][j][k + 1][2][3];
        lhs[i][j][k][2][2][4] = tmp2 * fjac[i][j][k + 1][2][4] - tmp1 * njac[i][j][k + 1][2][4];
        lhs[i][j][k][2][3][0] = tmp2 * fjac[i][j][k + 1][3][0] - tmp1 * njac[i][j][k + 1][3][0];
        lhs[i][j][k][2][3][1] = tmp2 * fjac[i][j][k + 1][3][1] - tmp1 * njac[i][j][k + 1][3][1];
        lhs[i][j][k][2][3][2] = tmp2 * fjac[i][j][k + 1][3][2] - tmp1 * njac[i][j][k + 1][3][2];
        lhs[i][j][k][2][3][3] = tmp2 * fjac[i][j][k + 1][3][3] - tmp1 * njac[i][j][k + 1][3][3] - tmp1 * dz4;
        lhs[i][j][k][2][3][4] = tmp2 * fjac[i][j][k + 1][3][4] - tmp1 * njac[i][j][k + 1][3][4];
        lhs[i][j][k][2][4][0] = tmp2 * fjac[i][j][k + 1][4][0] - tmp1 * njac[i][j][k + 1][4][0];
        lhs[i][j][k][2][4][1] = tmp2 * fjac[i][j][k + 1][4][1] - tmp1 * njac[i][j][k + 1][4][1];
        lhs[i][j][k][2][4][2] = tmp2 * fjac[i][j][k + 1][4][2] - tmp1 * njac[i][j][k + 1][4][2];
        lhs[i][j][k][2][4][3] = tmp2 * fjac[i][j][k + 1][4][3] - tmp1 * njac[i][j][k + 1][4][3];
        lhs[i][j][k][2][4][4] = tmp2 * fjac[i][j][k + 1][4][4] - tmp1 * njac[i][j][k + 1][4][4] - tmp1 * dz5;
      }
    }
  }
}

static void compute_rhs()
{
  int i;
  int j;
  int k;
  int m;
  double rho_inv;
  double uijk;
  double up1;
  double um1;
  double vijk;
  double vp1;
  double vm1;
  double wijk;
  double wp1;
  double wm1;
  
#pragma omp parallel for private (rho_inv,i,j,k)
  for (i = 0; i <= grid_points[0] - 1; i += 1) {
    
#pragma omp parallel for private (rho_inv,j,k)
    for (j = 0; j <= grid_points[1] - 1; j += 1) {
      
#pragma omp parallel for private (rho_inv,k)
      for (k = 0; k <= grid_points[2] - 1; k += 1) {
        rho_inv = 1.0 / u[i][j][k][0];
        rho_i[i][j][k] = rho_inv;
        us[i][j][k] = u[i][j][k][1] * rho_inv;
        vs[i][j][k] = u[i][j][k][2] * rho_inv;
        ws[i][j][k] = u[i][j][k][3] * rho_inv;
        square[i][j][k] = 0.5 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * rho_inv;
        qs[i][j][k] = square[i][j][k] * rho_inv;
      }
    }
  }
  
#pragma omp parallel for private (i,j,k,m)
  for (i = 0; i <= grid_points[0] - 1; i += 1) {
    
#pragma omp parallel for private (j,k,m)
    for (j = 0; j <= grid_points[1] - 1; j += 1) {
      
#pragma omp parallel for private (k,m)
      for (k = 0; k <= grid_points[2] - 1; k += 1) {
        
#pragma omp parallel for private (m)
        for (m = 0; m <= 4; m += 1) {
          rhs[i][j][k][m] = forcing[i][j][k][m];
        }
      }
    }
  }
  
#pragma omp parallel for private (uijk,up1,um1,i,j,k)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (uijk,up1,um1,j,k)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      
#pragma omp parallel for private (uijk,up1,um1,k) firstprivate (tx2,xxcon2,xxcon3,xxcon4,xxcon5,dx1tx1,dx2tx1,dx3tx1,dx4tx1,dx5tx1,c1,c2,con43)
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        uijk = us[i][j][k];
        up1 = us[i + 1][j][k];
        um1 = us[i - 1][j][k];
        rhs[i][j][k][0] = rhs[i][j][k][0] + dx1tx1 * (u[i + 1][j][k][0] - 2.0 * u[i][j][k][0] + u[i - 1][j][k][0]) - tx2 * (u[i + 1][j][k][1] - u[i - 1][j][k][1]);
        rhs[i][j][k][1] = rhs[i][j][k][1] + dx2tx1 * (u[i + 1][j][k][1] - 2.0 * u[i][j][k][1] + u[i - 1][j][k][1]) + xxcon2 * con43 * (up1 - 2.0 * uijk + um1) - tx2 * (u[i + 1][j][k][1] * up1 - u[i - 1][j][k][1] * um1 + (u[i + 1][j][k][4] - square[i + 1][j][k] - u[i - 1][j][k][4] + square[i - 1][j][k]) * c2);
        rhs[i][j][k][2] = rhs[i][j][k][2] + dx3tx1 * (u[i + 1][j][k][2] - 2.0 * u[i][j][k][2] + u[i - 1][j][k][2]) + xxcon2 * (vs[i + 1][j][k] - 2.0 * vs[i][j][k] + vs[i - 1][j][k]) - tx2 * (u[i + 1][j][k][2] * up1 - u[i - 1][j][k][2] * um1);
        rhs[i][j][k][3] = rhs[i][j][k][3] + dx4tx1 * (u[i + 1][j][k][3] - 2.0 * u[i][j][k][3] + u[i - 1][j][k][3]) + xxcon2 * (ws[i + 1][j][k] - 2.0 * ws[i][j][k] + ws[i - 1][j][k]) - tx2 * (u[i + 1][j][k][3] * up1 - u[i - 1][j][k][3] * um1);
        rhs[i][j][k][4] = rhs[i][j][k][4] + dx5tx1 * (u[i + 1][j][k][4] - 2.0 * u[i][j][k][4] + u[i - 1][j][k][4]) + xxcon3 * (qs[i + 1][j][k] - 2.0 * qs[i][j][k] + qs[i - 1][j][k]) + xxcon4 * (up1 * up1 - 2.0 * uijk * uijk + um1 * um1) + xxcon5 * (u[i + 1][j][k][4] * rho_i[i + 1][j][k] - 2.0 * u[i][j][k][4] * rho_i[i][j][k] + u[i - 1][j][k][4] * rho_i[i - 1][j][k]) - tx2 * ((c1 * u[i + 1][j][k][4] - c2 * square[i + 1][j][k]) * up1 - (c1 * u[i - 1][j][k][4] - c2 * square[i - 1][j][k]) * um1);
      }
    }
  }
  i = 1;
  
#pragma omp parallel for private (j,k,m)
  for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
    
#pragma omp parallel for private (k,m)
    for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
      
#pragma omp parallel for private (m) firstprivate (dssp,i)
      for (m = 0; m <= 4; m += 1) {
        rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (5.0 * u[i][j][k][m] - 4.0 * u[i + 1][j][k][m] + u[i + 2][j][k][m]);
      }
    }
  }
  i = 2;
  
#pragma omp parallel for private (j,k,m)
  for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
    
#pragma omp parallel for private (k,m)
    for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
      
#pragma omp parallel for private (m) firstprivate (dssp,i)
      for (m = 0; m <= 4; m += 1) {
        rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (- 4.0 * u[i - 1][j][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i + 1][j][k][m] + u[i + 2][j][k][m]);
      }
    }
  }
  
#pragma omp parallel for private (i,j,k,m)
  for (i = 3; i <= grid_points[0] - 3 - 1; i += 1) {
    
#pragma omp parallel for private (j,k,m)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      
#pragma omp parallel for private (k,m)
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        
#pragma omp parallel for private (m) firstprivate (dssp)
        for (m = 0; m <= 4; m += 1) {
          rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i - 2][j][k][m] - 4.0 * u[i - 1][j][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i + 1][j][k][m] + u[i + 2][j][k][m]);
        }
      }
    }
  }
  i = grid_points[0] - 3;
  
#pragma omp parallel for private (j,k,m)
  for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
    
#pragma omp parallel for private (k,m)
    for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
      
#pragma omp parallel for private (m) firstprivate (dssp,i)
      for (m = 0; m <= 4; m += 1) {
        rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i - 2][j][k][m] - 4.0 * u[i - 1][j][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i + 1][j][k][m]);
      }
    }
  }
  i = grid_points[0] - 2;
  
#pragma omp parallel for private (j,k,m)
  for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
    
#pragma omp parallel for private (k,m)
    for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
      
#pragma omp parallel for private (m) firstprivate (dssp,i)
      for (m = 0; m <= 4; m += 1) {
        rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i - 2][j][k][m] - 4. * u[i - 1][j][k][m] + 5.0 * u[i][j][k][m]);
      }
    }
  }
  
#pragma omp parallel for private (vijk,vp1,vm1,i,j,k)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (vijk,vp1,vm1,j,k)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      
#pragma omp parallel for private (vijk,vp1,vm1,k) firstprivate (ty2,yycon2,yycon3,yycon4,yycon5,dy1ty1,dy2ty1,dy3ty1,dy4ty1,dy5ty1,c1,c2,con43)
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        vijk = vs[i][j][k];
        vp1 = vs[i][j + 1][k];
        vm1 = vs[i][j - 1][k];
        rhs[i][j][k][0] = rhs[i][j][k][0] + dy1ty1 * (u[i][j + 1][k][0] - 2.0 * u[i][j][k][0] + u[i][j - 1][k][0]) - ty2 * (u[i][j + 1][k][2] - u[i][j - 1][k][2]);
        rhs[i][j][k][1] = rhs[i][j][k][1] + dy2ty1 * (u[i][j + 1][k][1] - 2.0 * u[i][j][k][1] + u[i][j - 1][k][1]) + yycon2 * (us[i][j + 1][k] - 2.0 * us[i][j][k] + us[i][j - 1][k]) - ty2 * (u[i][j + 1][k][1] * vp1 - u[i][j - 1][k][1] * vm1);
        rhs[i][j][k][2] = rhs[i][j][k][2] + dy3ty1 * (u[i][j + 1][k][2] - 2.0 * u[i][j][k][2] + u[i][j - 1][k][2]) + yycon2 * con43 * (vp1 - 2.0 * vijk + vm1) - ty2 * (u[i][j + 1][k][2] * vp1 - u[i][j - 1][k][2] * vm1 + (u[i][j + 1][k][4] - square[i][j + 1][k] - u[i][j - 1][k][4] + square[i][j - 1][k]) * c2);
        rhs[i][j][k][3] = rhs[i][j][k][3] + dy4ty1 * (u[i][j + 1][k][3] - 2.0 * u[i][j][k][3] + u[i][j - 1][k][3]) + yycon2 * (ws[i][j + 1][k] - 2.0 * ws[i][j][k] + ws[i][j - 1][k]) - ty2 * (u[i][j + 1][k][3] * vp1 - u[i][j - 1][k][3] * vm1);
        rhs[i][j][k][4] = rhs[i][j][k][4] + dy5ty1 * (u[i][j + 1][k][4] - 2.0 * u[i][j][k][4] + u[i][j - 1][k][4]) + yycon3 * (qs[i][j + 1][k] - 2.0 * qs[i][j][k] + qs[i][j - 1][k]) + yycon4 * (vp1 * vp1 - 2.0 * vijk * vijk + vm1 * vm1) + yycon5 * (u[i][j + 1][k][4] * rho_i[i][j + 1][k] - 2.0 * u[i][j][k][4] * rho_i[i][j][k] + u[i][j - 1][k][4] * rho_i[i][j - 1][k]) - ty2 * ((c1 * u[i][j + 1][k][4] - c2 * square[i][j + 1][k]) * vp1 - (c1 * u[i][j - 1][k][4] - c2 * square[i][j - 1][k]) * vm1);
      }
    }
  }
  j = 1;
  
#pragma omp parallel for private (i,k,m)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (k,m)
    for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
      
#pragma omp parallel for private (m) firstprivate (dssp,j)
      for (m = 0; m <= 4; m += 1) {
        rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (5.0 * u[i][j][k][m] - 4.0 * u[i][j + 1][k][m] + u[i][j + 2][k][m]);
      }
    }
  }
  j = 2;
  
#pragma omp parallel for private (i,k,m)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (k,m)
    for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
      
#pragma omp parallel for private (m) firstprivate (dssp,j)
      for (m = 0; m <= 4; m += 1) {
        rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (- 4.0 * u[i][j - 1][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j + 1][k][m] + u[i][j + 2][k][m]);
      }
    }
  }
  
#pragma omp parallel for private (i,j,k,m)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (j,k,m)
    for (j = 3; j <= grid_points[1] - 3 - 1; j += 1) {
      
#pragma omp parallel for private (k,m)
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        
#pragma omp parallel for private (m) firstprivate (dssp)
        for (m = 0; m <= 4; m += 1) {
          rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i][j - 2][k][m] - 4.0 * u[i][j - 1][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j + 1][k][m] + u[i][j + 2][k][m]);
        }
      }
    }
  }
  j = grid_points[1] - 3;
  
#pragma omp parallel for private (i,k,m)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (k,m)
    for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
      
#pragma omp parallel for private (m) firstprivate (dssp,j)
      for (m = 0; m <= 4; m += 1) {
        rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i][j - 2][k][m] - 4.0 * u[i][j - 1][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j + 1][k][m]);
      }
    }
  }
  j = grid_points[1] - 2;
  
#pragma omp parallel for private (i,k,m)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (k,m)
    for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
      
#pragma omp parallel for private (m) firstprivate (dssp,j)
      for (m = 0; m <= 4; m += 1) {
        rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i][j - 2][k][m] - 4. * u[i][j - 1][k][m] + 5. * u[i][j][k][m]);
      }
    }
  }
  
#pragma omp parallel for private (wijk,wp1,wm1,i,j,k)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (wijk,wp1,wm1,j,k)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      
#pragma omp parallel for private (wijk,wp1,wm1,k) firstprivate (tz2,zzcon2,zzcon3,zzcon4,zzcon5,dz1tz1,dz2tz1,dz3tz1,dz4tz1,dz5tz1,c1,c2,con43)
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        wijk = ws[i][j][k];
        wp1 = ws[i][j][k + 1];
        wm1 = ws[i][j][k - 1];
        rhs[i][j][k][0] = rhs[i][j][k][0] + dz1tz1 * (u[i][j][k + 1][0] - 2.0 * u[i][j][k][0] + u[i][j][k - 1][0]) - tz2 * (u[i][j][k + 1][3] - u[i][j][k - 1][3]);
        rhs[i][j][k][1] = rhs[i][j][k][1] + dz2tz1 * (u[i][j][k + 1][1] - 2.0 * u[i][j][k][1] + u[i][j][k - 1][1]) + zzcon2 * (us[i][j][k + 1] - 2.0 * us[i][j][k] + us[i][j][k - 1]) - tz2 * (u[i][j][k + 1][1] * wp1 - u[i][j][k - 1][1] * wm1);
        rhs[i][j][k][2] = rhs[i][j][k][2] + dz3tz1 * (u[i][j][k + 1][2] - 2.0 * u[i][j][k][2] + u[i][j][k - 1][2]) + zzcon2 * (vs[i][j][k + 1] - 2.0 * vs[i][j][k] + vs[i][j][k - 1]) - tz2 * (u[i][j][k + 1][2] * wp1 - u[i][j][k - 1][2] * wm1);
        rhs[i][j][k][3] = rhs[i][j][k][3] + dz4tz1 * (u[i][j][k + 1][3] - 2.0 * u[i][j][k][3] + u[i][j][k - 1][3]) + zzcon2 * con43 * (wp1 - 2.0 * wijk + wm1) - tz2 * (u[i][j][k + 1][3] * wp1 - u[i][j][k - 1][3] * wm1 + (u[i][j][k + 1][4] - square[i][j][k + 1] - u[i][j][k - 1][4] + square[i][j][k - 1]) * c2);
        rhs[i][j][k][4] = rhs[i][j][k][4] + dz5tz1 * (u[i][j][k + 1][4] - 2.0 * u[i][j][k][4] + u[i][j][k - 1][4]) + zzcon3 * (qs[i][j][k + 1] - 2.0 * qs[i][j][k] + qs[i][j][k - 1]) + zzcon4 * (wp1 * wp1 - 2.0 * wijk * wijk + wm1 * wm1) + zzcon5 * (u[i][j][k + 1][4] * rho_i[i][j][k + 1] - 2.0 * u[i][j][k][4] * rho_i[i][j][k] + u[i][j][k - 1][4] * rho_i[i][j][k - 1]) - tz2 * ((c1 * u[i][j][k + 1][4] - c2 * square[i][j][k + 1]) * wp1 - (c1 * u[i][j][k - 1][4] - c2 * square[i][j][k - 1]) * wm1);
      }
    }
  }
  k = 1;
  
#pragma omp parallel for private (i,j,m)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (j,m)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      
#pragma omp parallel for private (m) firstprivate (dssp,k)
      for (m = 0; m <= 4; m += 1) {
        rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (5.0 * u[i][j][k][m] - 4.0 * u[i][j][k + 1][m] + u[i][j][k + 2][m]);
      }
    }
  }
  k = 2;
  
#pragma omp parallel for private (i,j,m)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (j,m)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      
#pragma omp parallel for private (m) firstprivate (dssp,k)
      for (m = 0; m <= 4; m += 1) {
        rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (- 4.0 * u[i][j][k - 1][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j][k + 1][m] + u[i][j][k + 2][m]);
      }
    }
  }
  
#pragma omp parallel for private (i,j,k,m)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (j,k,m)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      
#pragma omp parallel for private (k,m)
      for (k = 3; k <= grid_points[2] - 3 - 1; k += 1) {
        
#pragma omp parallel for private (m) firstprivate (dssp)
        for (m = 0; m <= 4; m += 1) {
          rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i][j][k - 2][m] - 4.0 * u[i][j][k - 1][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j][k + 1][m] + u[i][j][k + 2][m]);
        }
      }
    }
  }
  k = grid_points[2] - 3;
  
#pragma omp parallel for private (i,j,m)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (j,m)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      
#pragma omp parallel for private (m) firstprivate (dssp,k)
      for (m = 0; m <= 4; m += 1) {
        rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i][j][k - 2][m] - 4.0 * u[i][j][k - 1][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j][k + 1][m]);
      }
    }
  }
  k = grid_points[2] - 2;
  
#pragma omp parallel for private (i,j,m)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (j,m)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      
#pragma omp parallel for private (m) firstprivate (dssp,k)
      for (m = 0; m <= 4; m += 1) {
        rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i][j][k - 2][m] - 4.0 * u[i][j][k - 1][m] + 5.0 * u[i][j][k][m]);
      }
    }
  }
  
#pragma omp parallel for private (i,j,k,m)
  for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
    
#pragma omp parallel for private (i,k,m)
    for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
      
#pragma omp parallel for private (i,m)
      for (m = 0; m <= 4; m += 1) {
        
#pragma omp parallel for private (i) firstprivate (dt)
        for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
          rhs[i][j][k][m] = rhs[i][j][k][m] * dt;
        }
      }
    }
  }
}

static void set_constants()
{
  ce[0][0] = 2.0;
  ce[0][1] = 0.0;
  ce[0][2] = 0.0;
  ce[0][3] = 4.0;
  ce[0][4] = 5.0;
  ce[0][5] = 3.0;
  ce[0][6] = 0.5;
  ce[0][7] = 0.02;
  ce[0][8] = 0.01;
  ce[0][9] = 0.03;
  ce[0][10] = 0.5;
  ce[0][11] = 0.4;
  ce[0][12] = 0.3;
  ce[1][0] = 1.0;
  ce[1][1] = 0.0;
  ce[1][2] = 0.0;
  ce[1][3] = 0.0;
  ce[1][4] = 1.0;
  ce[1][5] = 2.0;
  ce[1][6] = 3.0;
  ce[1][7] = 0.01;
  ce[1][8] = 0.03;
  ce[1][9] = 0.02;
  ce[1][10] = 0.4;
  ce[1][11] = 0.3;
  ce[1][12] = 0.5;
  ce[2][0] = 2.0;
  ce[2][1] = 2.0;
  ce[2][2] = 0.0;
  ce[2][3] = 0.0;
  ce[2][4] = 0.0;
  ce[2][5] = 2.0;
  ce[2][6] = 3.0;
  ce[2][7] = 0.04;
  ce[2][8] = 0.03;
  ce[2][9] = 0.05;
  ce[2][10] = 0.3;
  ce[2][11] = 0.5;
  ce[2][12] = 0.4;
  ce[3][0] = 2.0;
  ce[3][1] = 2.0;
  ce[3][2] = 0.0;
  ce[3][3] = 0.0;
  ce[3][4] = 0.0;
  ce[3][5] = 2.0;
  ce[3][6] = 3.0;
  ce[3][7] = 0.03;
  ce[3][8] = 0.05;
  ce[3][9] = 0.04;
  ce[3][10] = 0.2;
  ce[3][11] = 0.1;
  ce[3][12] = 0.3;
  ce[4][0] = 5.0;
  ce[4][1] = 4.0;
  ce[4][2] = 3.0;
  ce[4][3] = 2.0;
  ce[4][4] = 0.1;
  ce[4][5] = 0.4;
  ce[4][6] = 0.3;
  ce[4][7] = 0.05;
  ce[4][8] = 0.04;
  ce[4][9] = 0.03;
  ce[4][10] = 0.1;
  ce[4][11] = 0.3;
  ce[4][12] = 0.2;
  c1 = 1.4;
  c2 = 0.4;
  c3 = 0.1;
  c4 = 1.0;
  c5 = 1.4;
  dnxm1 = 1.0 / ((double )(grid_points[0] - 1));
  dnym1 = 1.0 / ((double )(grid_points[1] - 1));
  dnzm1 = 1.0 / ((double )(grid_points[2] - 1));
  c1c2 = c1 * c2;
  c1c5 = c1 * c5;
  c3c4 = c3 * c4;
  c1345 = c1c5 * c3c4;
  conz1 = 1.0 - c1c5;
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
  dxmax = (dx3 > dx4?dx3 : dx4);
  dymax = (dy2 > dy4?dy2 : dy4);
  dzmax = (dz2 > dz3?dz2 : dz3);
  dssp = 0.25 * ((dx1 > ((dy1 > dz1?dy1 : dz1))?dx1 : ((dy1 > dz1?dy1 : dz1))));
  c4dssp = 4.0 * dssp;
  c5dssp = 5.0 * dssp;
  dttx1 = dt * tx1;
  dttx2 = dt * tx2;
  dtty1 = dt * ty1;
  dtty2 = dt * ty2;
  dttz1 = dt * tz1;
  dttz2 = dt * tz2;
  c2dttx1 = 2.0 * dttx1;
  c2dtty1 = 2.0 * dtty1;
  c2dttz1 = 2.0 * dttz1;
  dtdssp = dt * dssp;
  comz1 = dtdssp;
  comz4 = 4.0 * dtdssp;
  comz5 = 5.0 * dtdssp;
  comz6 = 6.0 * dtdssp;
  c3c4tx3 = c3c4 * tx3;
  c3c4ty3 = c3c4 * ty3;
  c3c4tz3 = c3c4 * tz3;
  dx1tx1 = dx1 * tx1;
  dx2tx1 = dx2 * tx1;
  dx3tx1 = dx3 * tx1;
  dx4tx1 = dx4 * tx1;
  dx5tx1 = dx5 * tx1;
  dy1ty1 = dy1 * ty1;
  dy2ty1 = dy2 * ty1;
  dy3ty1 = dy3 * ty1;
  dy4ty1 = dy4 * ty1;
  dy5ty1 = dy5 * ty1;
  dz1tz1 = dz1 * tz1;
  dz2tz1 = dz2 * tz1;
  dz3tz1 = dz3 * tz1;
  dz4tz1 = dz4 * tz1;
  dz5tz1 = dz5 * tz1;
  c2iv = 2.5;
  con43 = 4.0 / 3.0;
  con16 = 1.0 / 6.0;
  xxcon1 = c3c4tx3 * con43 * tx3;
  xxcon2 = c3c4tx3 * tx3;
  xxcon3 = c3c4tx3 * conz1 * tx3;
  xxcon4 = c3c4tx3 * con16 * tx3;
  xxcon5 = c3c4tx3 * c1c5 * tx3;
  yycon1 = c3c4ty3 * con43 * ty3;
  yycon2 = c3c4ty3 * ty3;
  yycon3 = c3c4ty3 * conz1 * ty3;
  yycon4 = c3c4ty3 * con16 * ty3;
  yycon5 = c3c4ty3 * c1c5 * ty3;
  zzcon1 = c3c4tz3 * con43 * tz3;
  zzcon2 = c3c4tz3 * tz3;
  zzcon3 = c3c4tz3 * conz1 * tz3;
  zzcon4 = c3c4tz3 * con16 * tz3;
  zzcon5 = c3c4tz3 * c1c5 * tz3;
}

static void verify(int no_time_steps,char *class,boolean *verified)
{
  double xcrref[5];
  double xceref[5];
  double xcrdif[5];
  double xcedif[5];
  double epsilon;
  double xce[5];
  double xcr[5];
  double dtref;
  int m;
  epsilon = 1.0e-08;
  error_norm(xce);
  compute_rhs();
  rhs_norm(xcr);
  
#pragma omp parallel for private (m)
  for (m = 0; m <= 4; m += 1) {
    xcr[m] = xcr[m] / dt;
  }
   *class = 'U';
   *verified = 1;
  
#pragma omp parallel for private (m)
  for (m = 0; m <= 4; m += 1) {
    xcrref[m] = 1.0;
    xceref[m] = 1.0;
  }
  if (grid_points[0] == 12 && grid_points[1] == 12 && grid_points[2] == 12 && no_time_steps == 60) {
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
  }
   else if (grid_points[0] == 24 && grid_points[1] == 24 && grid_points[2] == 24 && no_time_steps == 200) {
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
  }
   else if (grid_points[0] == 64 && grid_points[1] == 64 && grid_points[2] == 64 && no_time_steps == 200) {
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
  }
   else if (grid_points[0] == 102 && grid_points[1] == 102 && grid_points[2] == 102 && no_time_steps == 200) {
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
  }
   else if (grid_points[0] == 162 && grid_points[1] == 162 && grid_points[2] == 162 && no_time_steps == 200) {
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
  }
   else {
     *verified = 0;
  }
  for (m = 0; m <= 4; m += 1) {
    xcrdif[m] = fabs((xcr[m] - xcrref[m]) / xcrref[m]);
    xcedif[m] = fabs((xce[m] - xceref[m]) / xceref[m]);
  }
  if (( *class) != 'U') {
    printf(" Verification being performed for class %1c\n",( *class));
    printf(" accuracy setting for epsilon = %20.13e\n",epsilon);
    if (fabs(dt - dtref) > epsilon) {
       *verified = 0;
       *class = 'U';
      printf(" DT does not match the reference value of %15.8e\n",dtref);
    }
  }
   else {
    printf(" Unknown class\n");
  }
  if (( *class) != 'U') {
    printf(" Comparison of RMS-norms of residual\n");
  }
   else {
    printf(" RMS-norms of residual\n");
  }
  for (m = 0; m <= 4; m += 1) {
    if (( *class) == 'U') {
      printf("          %2d%20.13e\n",m,xcr[m]);
    }
     else if (xcrdif[m] > epsilon) {
       *verified = 0;
      printf(" FAILURE: %2d%20.13e%20.13e%20.13e\n",m,xcr[m],xcrref[m],xcrdif[m]);
    }
     else {
      printf("          %2d%20.13e%20.13e%20.13e\n",m,xcr[m],xcrref[m],xcrdif[m]);
    }
  }
  if (( *class) != 'U') {
    printf(" Comparison of RMS-norms of solution error\n");
  }
   else {
    printf(" RMS-norms of solution error\n");
  }
  for (m = 0; m <= 4; m += 1) {
    if (( *class) == 'U') {
      printf("          %2d%20.13e\n",m,xce[m]);
    }
     else if (xcedif[m] > epsilon) {
       *verified = 0;
      printf(" FAILURE: %2d%20.13e%20.13e%20.13e\n",m,xce[m],xceref[m],xcedif[m]);
    }
     else {
      printf("          %2d%20.13e%20.13e%20.13e\n",m,xce[m],xceref[m],xcedif[m]);
    }
  }
  if (( *class) == 'U') {
    printf(" No reference values provided\n");
    printf(" No verification performed\n");
  }
   else if ( *verified == 1) {
    printf(" Verification Successful\n");
  }
   else {
    printf(" Verification failed\n");
  }
}

static void x_solve()
{
  lhsx();
  x_solve_cell();
  x_backsubstitute();
}

static void x_backsubstitute()
{
  int i;
  int j;
  int k;
  int m;
  int n;
  for (i = grid_points[0] - 2; i >= 0; i += -1) {
    
#pragma omp parallel for private (j,k,m,n)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      
#pragma omp parallel for private (k,m,n)
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        
#pragma omp parallel for private (m,n)
        for (m = 0; m <= 4; m += 1) {
          for (n = 0; n <= 4; n += 1) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - lhs[i][j][k][2][m][n] * rhs[i + 1][j][k][n];
          }
        }
      }
    }
  }
}

static void x_solve_cell()
{
  int i;
  int j;
  int k;
  int isize;
  isize = grid_points[0] - 1;
  for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
    for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
      binvcrhs(lhs[0][j][k][1],lhs[0][j][k][2],rhs[0][j][k]);
    }
  }
  for (i = 1; i <= isize - 1; i += 1) {
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        matvec_sub(lhs[i][j][k][0],rhs[i - 1][j][k],rhs[i][j][k]);
        matmul_sub(lhs[i][j][k][0],lhs[i - 1][j][k][2],lhs[i][j][k][1]);
        binvcrhs(lhs[i][j][k][1],lhs[i][j][k][2],rhs[i][j][k]);
      }
    }
  }
  for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
    for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
      matvec_sub(lhs[isize][j][k][0],rhs[isize - 1][j][k],rhs[isize][j][k]);
      matmul_sub(lhs[isize][j][k][0],lhs[isize - 1][j][k][2],lhs[isize][j][k][1]);
      binvrhs(lhs[i][j][k][1],rhs[i][j][k]);
    }
  }
}

static void matvec_sub(double ablock[5][5],double avec[5],double bvec[5])
{
  int i;
  for (i = 0; i <= 4; i += 1) {
    bvec[i] = bvec[i] - ablock[i][0] * avec[0] - ablock[i][1] * avec[1] - ablock[i][2] * avec[2] - ablock[i][3] * avec[3] - ablock[i][4] * avec[4];
  }
}

static void matmul_sub(double ablock[5][5],double bblock[5][5],double cblock[5][5])
{
  int j;
  for (j = 0; j <= 4; j += 1) {
    cblock[0][j] = cblock[0][j] - ablock[0][0] * bblock[0][j] - ablock[0][1] * bblock[1][j] - ablock[0][2] * bblock[2][j] - ablock[0][3] * bblock[3][j] - ablock[0][4] * bblock[4][j];
    cblock[1][j] = cblock[1][j] - ablock[1][0] * bblock[0][j] - ablock[1][1] * bblock[1][j] - ablock[1][2] * bblock[2][j] - ablock[1][3] * bblock[3][j] - ablock[1][4] * bblock[4][j];
    cblock[2][j] = cblock[2][j] - ablock[2][0] * bblock[0][j] - ablock[2][1] * bblock[1][j] - ablock[2][2] * bblock[2][j] - ablock[2][3] * bblock[3][j] - ablock[2][4] * bblock[4][j];
    cblock[3][j] = cblock[3][j] - ablock[3][0] * bblock[0][j] - ablock[3][1] * bblock[1][j] - ablock[3][2] * bblock[2][j] - ablock[3][3] * bblock[3][j] - ablock[3][4] * bblock[4][j];
    cblock[4][j] = cblock[4][j] - ablock[4][0] * bblock[0][j] - ablock[4][1] * bblock[1][j] - ablock[4][2] * bblock[2][j] - ablock[4][3] * bblock[3][j] - ablock[4][4] * bblock[4][j];
  }
}

static void binvcrhs(double lhs[5][5],double c[5][5],double r[5])
{
  double pivot;
  double coeff;
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
  coeff = lhs[1][0];
  lhs[1][1] = lhs[1][1] - coeff * lhs[0][1];
  lhs[1][2] = lhs[1][2] - coeff * lhs[0][2];
  lhs[1][3] = lhs[1][3] - coeff * lhs[0][3];
  lhs[1][4] = lhs[1][4] - coeff * lhs[0][4];
  c[1][0] = c[1][0] - coeff * c[0][0];
  c[1][1] = c[1][1] - coeff * c[0][1];
  c[1][2] = c[1][2] - coeff * c[0][2];
  c[1][3] = c[1][3] - coeff * c[0][3];
  c[1][4] = c[1][4] - coeff * c[0][4];
  r[1] = r[1] - coeff * r[0];
  coeff = lhs[2][0];
  lhs[2][1] = lhs[2][1] - coeff * lhs[0][1];
  lhs[2][2] = lhs[2][2] - coeff * lhs[0][2];
  lhs[2][3] = lhs[2][3] - coeff * lhs[0][3];
  lhs[2][4] = lhs[2][4] - coeff * lhs[0][4];
  c[2][0] = c[2][0] - coeff * c[0][0];
  c[2][1] = c[2][1] - coeff * c[0][1];
  c[2][2] = c[2][2] - coeff * c[0][2];
  c[2][3] = c[2][3] - coeff * c[0][3];
  c[2][4] = c[2][4] - coeff * c[0][4];
  r[2] = r[2] - coeff * r[0];
  coeff = lhs[3][0];
  lhs[3][1] = lhs[3][1] - coeff * lhs[0][1];
  lhs[3][2] = lhs[3][2] - coeff * lhs[0][2];
  lhs[3][3] = lhs[3][3] - coeff * lhs[0][3];
  lhs[3][4] = lhs[3][4] - coeff * lhs[0][4];
  c[3][0] = c[3][0] - coeff * c[0][0];
  c[3][1] = c[3][1] - coeff * c[0][1];
  c[3][2] = c[3][2] - coeff * c[0][2];
  c[3][3] = c[3][3] - coeff * c[0][3];
  c[3][4] = c[3][4] - coeff * c[0][4];
  r[3] = r[3] - coeff * r[0];
  coeff = lhs[4][0];
  lhs[4][1] = lhs[4][1] - coeff * lhs[0][1];
  lhs[4][2] = lhs[4][2] - coeff * lhs[0][2];
  lhs[4][3] = lhs[4][3] - coeff * lhs[0][3];
  lhs[4][4] = lhs[4][4] - coeff * lhs[0][4];
  c[4][0] = c[4][0] - coeff * c[0][0];
  c[4][1] = c[4][1] - coeff * c[0][1];
  c[4][2] = c[4][2] - coeff * c[0][2];
  c[4][3] = c[4][3] - coeff * c[0][3];
  c[4][4] = c[4][4] - coeff * c[0][4];
  r[4] = r[4] - coeff * r[0];
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
  coeff = lhs[0][1];
  lhs[0][2] = lhs[0][2] - coeff * lhs[1][2];
  lhs[0][3] = lhs[0][3] - coeff * lhs[1][3];
  lhs[0][4] = lhs[0][4] - coeff * lhs[1][4];
  c[0][0] = c[0][0] - coeff * c[1][0];
  c[0][1] = c[0][1] - coeff * c[1][1];
  c[0][2] = c[0][2] - coeff * c[1][2];
  c[0][3] = c[0][3] - coeff * c[1][3];
  c[0][4] = c[0][4] - coeff * c[1][4];
  r[0] = r[0] - coeff * r[1];
  coeff = lhs[2][1];
  lhs[2][2] = lhs[2][2] - coeff * lhs[1][2];
  lhs[2][3] = lhs[2][3] - coeff * lhs[1][3];
  lhs[2][4] = lhs[2][4] - coeff * lhs[1][4];
  c[2][0] = c[2][0] - coeff * c[1][0];
  c[2][1] = c[2][1] - coeff * c[1][1];
  c[2][2] = c[2][2] - coeff * c[1][2];
  c[2][3] = c[2][3] - coeff * c[1][3];
  c[2][4] = c[2][4] - coeff * c[1][4];
  r[2] = r[2] - coeff * r[1];
  coeff = lhs[3][1];
  lhs[3][2] = lhs[3][2] - coeff * lhs[1][2];
  lhs[3][3] = lhs[3][3] - coeff * lhs[1][3];
  lhs[3][4] = lhs[3][4] - coeff * lhs[1][4];
  c[3][0] = c[3][0] - coeff * c[1][0];
  c[3][1] = c[3][1] - coeff * c[1][1];
  c[3][2] = c[3][2] - coeff * c[1][2];
  c[3][3] = c[3][3] - coeff * c[1][3];
  c[3][4] = c[3][4] - coeff * c[1][4];
  r[3] = r[3] - coeff * r[1];
  coeff = lhs[4][1];
  lhs[4][2] = lhs[4][2] - coeff * lhs[1][2];
  lhs[4][3] = lhs[4][3] - coeff * lhs[1][3];
  lhs[4][4] = lhs[4][4] - coeff * lhs[1][4];
  c[4][0] = c[4][0] - coeff * c[1][0];
  c[4][1] = c[4][1] - coeff * c[1][1];
  c[4][2] = c[4][2] - coeff * c[1][2];
  c[4][3] = c[4][3] - coeff * c[1][3];
  c[4][4] = c[4][4] - coeff * c[1][4];
  r[4] = r[4] - coeff * r[1];
  pivot = 1.00 / lhs[2][2];
  lhs[2][3] = lhs[2][3] * pivot;
  lhs[2][4] = lhs[2][4] * pivot;
  c[2][0] = c[2][0] * pivot;
  c[2][1] = c[2][1] * pivot;
  c[2][2] = c[2][2] * pivot;
  c[2][3] = c[2][3] * pivot;
  c[2][4] = c[2][4] * pivot;
  r[2] = r[2] * pivot;
  coeff = lhs[0][2];
  lhs[0][3] = lhs[0][3] - coeff * lhs[2][3];
  lhs[0][4] = lhs[0][4] - coeff * lhs[2][4];
  c[0][0] = c[0][0] - coeff * c[2][0];
  c[0][1] = c[0][1] - coeff * c[2][1];
  c[0][2] = c[0][2] - coeff * c[2][2];
  c[0][3] = c[0][3] - coeff * c[2][3];
  c[0][4] = c[0][4] - coeff * c[2][4];
  r[0] = r[0] - coeff * r[2];
  coeff = lhs[1][2];
  lhs[1][3] = lhs[1][3] - coeff * lhs[2][3];
  lhs[1][4] = lhs[1][4] - coeff * lhs[2][4];
  c[1][0] = c[1][0] - coeff * c[2][0];
  c[1][1] = c[1][1] - coeff * c[2][1];
  c[1][2] = c[1][2] - coeff * c[2][2];
  c[1][3] = c[1][3] - coeff * c[2][3];
  c[1][4] = c[1][4] - coeff * c[2][4];
  r[1] = r[1] - coeff * r[2];
  coeff = lhs[3][2];
  lhs[3][3] = lhs[3][3] - coeff * lhs[2][3];
  lhs[3][4] = lhs[3][4] - coeff * lhs[2][4];
  c[3][0] = c[3][0] - coeff * c[2][0];
  c[3][1] = c[3][1] - coeff * c[2][1];
  c[3][2] = c[3][2] - coeff * c[2][2];
  c[3][3] = c[3][3] - coeff * c[2][3];
  c[3][4] = c[3][4] - coeff * c[2][4];
  r[3] = r[3] - coeff * r[2];
  coeff = lhs[4][2];
  lhs[4][3] = lhs[4][3] - coeff * lhs[2][3];
  lhs[4][4] = lhs[4][4] - coeff * lhs[2][4];
  c[4][0] = c[4][0] - coeff * c[2][0];
  c[4][1] = c[4][1] - coeff * c[2][1];
  c[4][2] = c[4][2] - coeff * c[2][2];
  c[4][3] = c[4][3] - coeff * c[2][3];
  c[4][4] = c[4][4] - coeff * c[2][4];
  r[4] = r[4] - coeff * r[2];
  pivot = 1.00 / lhs[3][3];
  lhs[3][4] = lhs[3][4] * pivot;
  c[3][0] = c[3][0] * pivot;
  c[3][1] = c[3][1] * pivot;
  c[3][2] = c[3][2] * pivot;
  c[3][3] = c[3][3] * pivot;
  c[3][4] = c[3][4] * pivot;
  r[3] = r[3] * pivot;
  coeff = lhs[0][3];
  lhs[0][4] = lhs[0][4] - coeff * lhs[3][4];
  c[0][0] = c[0][0] - coeff * c[3][0];
  c[0][1] = c[0][1] - coeff * c[3][1];
  c[0][2] = c[0][2] - coeff * c[3][2];
  c[0][3] = c[0][3] - coeff * c[3][3];
  c[0][4] = c[0][4] - coeff * c[3][4];
  r[0] = r[0] - coeff * r[3];
  coeff = lhs[1][3];
  lhs[1][4] = lhs[1][4] - coeff * lhs[3][4];
  c[1][0] = c[1][0] - coeff * c[3][0];
  c[1][1] = c[1][1] - coeff * c[3][1];
  c[1][2] = c[1][2] - coeff * c[3][2];
  c[1][3] = c[1][3] - coeff * c[3][3];
  c[1][4] = c[1][4] - coeff * c[3][4];
  r[1] = r[1] - coeff * r[3];
  coeff = lhs[2][3];
  lhs[2][4] = lhs[2][4] - coeff * lhs[3][4];
  c[2][0] = c[2][0] - coeff * c[3][0];
  c[2][1] = c[2][1] - coeff * c[3][1];
  c[2][2] = c[2][2] - coeff * c[3][2];
  c[2][3] = c[2][3] - coeff * c[3][3];
  c[2][4] = c[2][4] - coeff * c[3][4];
  r[2] = r[2] - coeff * r[3];
  coeff = lhs[4][3];
  lhs[4][4] = lhs[4][4] - coeff * lhs[3][4];
  c[4][0] = c[4][0] - coeff * c[3][0];
  c[4][1] = c[4][1] - coeff * c[3][1];
  c[4][2] = c[4][2] - coeff * c[3][2];
  c[4][3] = c[4][3] - coeff * c[3][3];
  c[4][4] = c[4][4] - coeff * c[3][4];
  r[4] = r[4] - coeff * r[3];
  pivot = 1.00 / lhs[4][4];
  c[4][0] = c[4][0] * pivot;
  c[4][1] = c[4][1] * pivot;
  c[4][2] = c[4][2] * pivot;
  c[4][3] = c[4][3] * pivot;
  c[4][4] = c[4][4] * pivot;
  r[4] = r[4] * pivot;
  coeff = lhs[0][4];
  c[0][0] = c[0][0] - coeff * c[4][0];
  c[0][1] = c[0][1] - coeff * c[4][1];
  c[0][2] = c[0][2] - coeff * c[4][2];
  c[0][3] = c[0][3] - coeff * c[4][3];
  c[0][4] = c[0][4] - coeff * c[4][4];
  r[0] = r[0] - coeff * r[4];
  coeff = lhs[1][4];
  c[1][0] = c[1][0] - coeff * c[4][0];
  c[1][1] = c[1][1] - coeff * c[4][1];
  c[1][2] = c[1][2] - coeff * c[4][2];
  c[1][3] = c[1][3] - coeff * c[4][3];
  c[1][4] = c[1][4] - coeff * c[4][4];
  r[1] = r[1] - coeff * r[4];
  coeff = lhs[2][4];
  c[2][0] = c[2][0] - coeff * c[4][0];
  c[2][1] = c[2][1] - coeff * c[4][1];
  c[2][2] = c[2][2] - coeff * c[4][2];
  c[2][3] = c[2][3] - coeff * c[4][3];
  c[2][4] = c[2][4] - coeff * c[4][4];
  r[2] = r[2] - coeff * r[4];
  coeff = lhs[3][4];
  c[3][0] = c[3][0] - coeff * c[4][0];
  c[3][1] = c[3][1] - coeff * c[4][1];
  c[3][2] = c[3][2] - coeff * c[4][2];
  c[3][3] = c[3][3] - coeff * c[4][3];
  c[3][4] = c[3][4] - coeff * c[4][4];
  r[3] = r[3] - coeff * r[4];
}

static void binvrhs(double lhs[5][5],double r[5])
{
  double pivot;
  double coeff;
  pivot = 1.00 / lhs[0][0];
  lhs[0][1] = lhs[0][1] * pivot;
  lhs[0][2] = lhs[0][2] * pivot;
  lhs[0][3] = lhs[0][3] * pivot;
  lhs[0][4] = lhs[0][4] * pivot;
  r[0] = r[0] * pivot;
  coeff = lhs[1][0];
  lhs[1][1] = lhs[1][1] - coeff * lhs[0][1];
  lhs[1][2] = lhs[1][2] - coeff * lhs[0][2];
  lhs[1][3] = lhs[1][3] - coeff * lhs[0][3];
  lhs[1][4] = lhs[1][4] - coeff * lhs[0][4];
  r[1] = r[1] - coeff * r[0];
  coeff = lhs[2][0];
  lhs[2][1] = lhs[2][1] - coeff * lhs[0][1];
  lhs[2][2] = lhs[2][2] - coeff * lhs[0][2];
  lhs[2][3] = lhs[2][3] - coeff * lhs[0][3];
  lhs[2][4] = lhs[2][4] - coeff * lhs[0][4];
  r[2] = r[2] - coeff * r[0];
  coeff = lhs[3][0];
  lhs[3][1] = lhs[3][1] - coeff * lhs[0][1];
  lhs[3][2] = lhs[3][2] - coeff * lhs[0][2];
  lhs[3][3] = lhs[3][3] - coeff * lhs[0][3];
  lhs[3][4] = lhs[3][4] - coeff * lhs[0][4];
  r[3] = r[3] - coeff * r[0];
  coeff = lhs[4][0];
  lhs[4][1] = lhs[4][1] - coeff * lhs[0][1];
  lhs[4][2] = lhs[4][2] - coeff * lhs[0][2];
  lhs[4][3] = lhs[4][3] - coeff * lhs[0][3];
  lhs[4][4] = lhs[4][4] - coeff * lhs[0][4];
  r[4] = r[4] - coeff * r[0];
  pivot = 1.00 / lhs[1][1];
  lhs[1][2] = lhs[1][2] * pivot;
  lhs[1][3] = lhs[1][3] * pivot;
  lhs[1][4] = lhs[1][4] * pivot;
  r[1] = r[1] * pivot;
  coeff = lhs[0][1];
  lhs[0][2] = lhs[0][2] - coeff * lhs[1][2];
  lhs[0][3] = lhs[0][3] - coeff * lhs[1][3];
  lhs[0][4] = lhs[0][4] - coeff * lhs[1][4];
  r[0] = r[0] - coeff * r[1];
  coeff = lhs[2][1];
  lhs[2][2] = lhs[2][2] - coeff * lhs[1][2];
  lhs[2][3] = lhs[2][3] - coeff * lhs[1][3];
  lhs[2][4] = lhs[2][4] - coeff * lhs[1][4];
  r[2] = r[2] - coeff * r[1];
  coeff = lhs[3][1];
  lhs[3][2] = lhs[3][2] - coeff * lhs[1][2];
  lhs[3][3] = lhs[3][3] - coeff * lhs[1][3];
  lhs[3][4] = lhs[3][4] - coeff * lhs[1][4];
  r[3] = r[3] - coeff * r[1];
  coeff = lhs[4][1];
  lhs[4][2] = lhs[4][2] - coeff * lhs[1][2];
  lhs[4][3] = lhs[4][3] - coeff * lhs[1][3];
  lhs[4][4] = lhs[4][4] - coeff * lhs[1][4];
  r[4] = r[4] - coeff * r[1];
  pivot = 1.00 / lhs[2][2];
  lhs[2][3] = lhs[2][3] * pivot;
  lhs[2][4] = lhs[2][4] * pivot;
  r[2] = r[2] * pivot;
  coeff = lhs[0][2];
  lhs[0][3] = lhs[0][3] - coeff * lhs[2][3];
  lhs[0][4] = lhs[0][4] - coeff * lhs[2][4];
  r[0] = r[0] - coeff * r[2];
  coeff = lhs[1][2];
  lhs[1][3] = lhs[1][3] - coeff * lhs[2][3];
  lhs[1][4] = lhs[1][4] - coeff * lhs[2][4];
  r[1] = r[1] - coeff * r[2];
  coeff = lhs[3][2];
  lhs[3][3] = lhs[3][3] - coeff * lhs[2][3];
  lhs[3][4] = lhs[3][4] - coeff * lhs[2][4];
  r[3] = r[3] - coeff * r[2];
  coeff = lhs[4][2];
  lhs[4][3] = lhs[4][3] - coeff * lhs[2][3];
  lhs[4][4] = lhs[4][4] - coeff * lhs[2][4];
  r[4] = r[4] - coeff * r[2];
  pivot = 1.00 / lhs[3][3];
  lhs[3][4] = lhs[3][4] * pivot;
  r[3] = r[3] * pivot;
  coeff = lhs[0][3];
  lhs[0][4] = lhs[0][4] - coeff * lhs[3][4];
  r[0] = r[0] - coeff * r[3];
  coeff = lhs[1][3];
  lhs[1][4] = lhs[1][4] - coeff * lhs[3][4];
  r[1] = r[1] - coeff * r[3];
  coeff = lhs[2][3];
  lhs[2][4] = lhs[2][4] - coeff * lhs[3][4];
  r[2] = r[2] - coeff * r[3];
  coeff = lhs[4][3];
  lhs[4][4] = lhs[4][4] - coeff * lhs[3][4];
  r[4] = r[4] - coeff * r[3];
  pivot = 1.00 / lhs[4][4];
  r[4] = r[4] * pivot;
  coeff = lhs[0][4];
  r[0] = r[0] - coeff * r[4];
  coeff = lhs[1][4];
  r[1] = r[1] - coeff * r[4];
  coeff = lhs[2][4];
  r[2] = r[2] - coeff * r[4];
  coeff = lhs[3][4];
  r[3] = r[3] - coeff * r[4];
}

static void y_solve()
{
  lhsy();
  y_solve_cell();
  y_backsubstitute();
}

static void y_backsubstitute()
{
  int i;
  int j;
  int k;
  int m;
  int n;
  for (j = grid_points[1] - 2; j >= 0; j += -1) {
    
#pragma omp parallel for private (i,k,m,n)
    for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
      
#pragma omp parallel for private (k,m,n)
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        
#pragma omp parallel for private (m,n)
        for (m = 0; m <= 4; m += 1) {
          for (n = 0; n <= 4; n += 1) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - lhs[i][j][k][2][m][n] * rhs[i][j + 1][k][n];
          }
        }
      }
    }
  }
}

static void y_solve_cell()
{
  int i;
  int j;
  int k;
  int jsize;
  jsize = grid_points[1] - 1;
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
      binvcrhs(lhs[i][0][k][1],lhs[i][0][k][2],rhs[i][0][k]);
    }
  }
  for (j = 1; j <= jsize - 1; j += 1) {
    for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
      for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
        matvec_sub(lhs[i][j][k][0],rhs[i][j - 1][k],rhs[i][j][k]);
        matmul_sub(lhs[i][j][k][0],lhs[i][j - 1][k][2],lhs[i][j][k][1]);
        binvcrhs(lhs[i][j][k][1],lhs[i][j][k][2],rhs[i][j][k]);
      }
    }
  }
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    for (k = 1; k <= grid_points[2] - 1 - 1; k += 1) {
      matvec_sub(lhs[i][jsize][k][0],rhs[i][jsize - 1][k],rhs[i][jsize][k]);
      matmul_sub(lhs[i][jsize][k][0],lhs[i][jsize - 1][k][2],lhs[i][jsize][k][1]);
      binvrhs(lhs[i][jsize][k][1],rhs[i][jsize][k]);
    }
  }
}

static void z_solve()
{
  lhsz();
  z_solve_cell();
  z_backsubstitute();
}

static void z_backsubstitute()
{
  int i;
  int j;
  int k;
  int m;
  int n;
  
#pragma omp parallel for private (i,j,k,m,n)
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    
#pragma omp parallel for private (j,k,m,n)
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      for (k = grid_points[2] - 2; k >= 0; k += -1) {
        
#pragma omp parallel for private (m,n)
        for (m = 0; m <= 4; m += 1) {
          for (n = 0; n <= 4; n += 1) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - lhs[i][j][k][2][m][n] * rhs[i][j][k + 1][n];
          }
        }
      }
    }
  }
}

static void z_solve_cell()
{
  int i;
  int j;
  int k;
  int ksize;
  ksize = grid_points[2] - 1;
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      binvcrhs(lhs[i][j][0][1],lhs[i][j][0][2],rhs[i][j][0]);
    }
  }
  for (k = 1; k <= ksize - 1; k += 1) {
    for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
      for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
        matvec_sub(lhs[i][j][k][0],rhs[i][j][k - 1],rhs[i][j][k]);
        matmul_sub(lhs[i][j][k][0],lhs[i][j][k - 1][2],lhs[i][j][k][1]);
        binvcrhs(lhs[i][j][k][1],lhs[i][j][k][2],rhs[i][j][k]);
      }
    }
  }
  for (i = 1; i <= grid_points[0] - 1 - 1; i += 1) {
    for (j = 1; j <= grid_points[1] - 1 - 1; j += 1) {
      matvec_sub(lhs[i][j][ksize][0],rhs[i][j][ksize - 1],rhs[i][j][ksize]);
      matmul_sub(lhs[i][j][ksize][0],lhs[i][j][ksize - 1][2],lhs[i][j][ksize][1]);
      binvrhs(lhs[i][j][ksize][1],rhs[i][j][ksize]);
    }
  }
}
