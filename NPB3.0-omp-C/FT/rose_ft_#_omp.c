#include "npb-C.h"
#include "global.h"
#include <omp.h> 
static void evolve(dcomplex u0[32][128][128],dcomplex u1[32][128][128],int t,int indexmap[32][128][128],int d[3]);
static void compute_initial_conditions(dcomplex u0[32][128][128],int d[3]);
static void ipow46(double a,int exponent,double *result);
static void setup();
static void compute_indexmap(int indexmap[32][128][128],int d[3]);
static void print_timers();
static void fft(int dir,dcomplex x1[32][128][128],dcomplex x2[32][128][128]);
static void cffts1(int is,int d[3],dcomplex x[32][128][128],dcomplex xout[32][128][128],dcomplex y0[128][18],dcomplex y1[128][18]);
static void cffts2(int is,int d[3],dcomplex x[32][128][128],dcomplex xout[32][128][128],dcomplex y0[128][18],dcomplex y1[128][18]);
static void cffts3(int is,int d[3],dcomplex x[32][128][128],dcomplex xout[32][128][128],dcomplex y0[128][18],dcomplex y1[128][18]);
static void fft_init(int n);
static void cfftz(int is,int m,int n,dcomplex x[128][18],dcomplex y[128][18]);
static void fftz2(int is,int l,int m,int n,int ny,int ny1,dcomplex u[128],dcomplex x[128][18],dcomplex y[128][18]);
static int ilog2(int n);
static void checksum(int i,dcomplex u1[32][128][128],int d[3]);
static void verify(int d1,int d2,int d3,int nt,boolean *verified,char *class);

int main(int argc,char **argv)
{
  int i;
  int ierr;
  static dcomplex u0[32][128][128];
  static dcomplex pad1[3];
  static dcomplex u1[32][128][128];
  static dcomplex pad2[3];
  static dcomplex u2[32][128][128];
  static dcomplex pad3[3];
  static int indexmap[32][128][128];
  int iter;
  int nthreads = 1;
  double total_time;
  double mflops;
  boolean verified;
  char class;
  for (i = 0; i <= 6; i += 1) {
    timer_clear(i);
  }
  setup();
  compute_indexmap(indexmap,dims[2]);
  compute_initial_conditions(u1,dims[0]);
  fft_init(dims[0][0]);
  fft(1,u1,u0);
  for (i = 0; i <= 6; i += 1) {
    timer_clear(i);
  }
  timer_start(0);
  if (0 == 1) 
    timer_start(1);
  compute_indexmap(indexmap,dims[2]);
  compute_initial_conditions(u1,dims[0]);
  fft_init(dims[0][0]);
  if (0 == 1) {
    timer_stop(1);
  }
  if (0 == 1) {
    timer_start(2);
  }
  fft(1,u1,u0);
  if (0 == 1) {
    timer_stop(2);
  }
  for (iter = 1; iter <= niter; iter += 1) {
    if (0 == 1) {
      timer_start(3);
    }
    evolve(u0,u1,iter,indexmap,dims[0]);
    if (0 == 1) {
      timer_stop(3);
    }
    if (0 == 1) {
      timer_start(2);
    }
    fft(- 1,u1,u2);
    if (0 == 1) {
      timer_stop(2);
    }
    if (0 == 1) {
      timer_start(4);
    }
    checksum(iter,u2,dims[0]);
    if (0 == 1) {
      timer_stop(4);
    }
  }
  verify(128,128,32,niter,&verified,&class);
{
#if defined(_OPENMP)
#endif     
  }
  timer_stop(0);
  total_time = timer_read(0);
  if (total_time != 0.0) {
    mflops = 1.0e-6 * ((double )524288) * (14.8157 + 7.19641 * log((double )524288) + (5.23518 + 7.21113 * log((double )524288)) * niter) / total_time;
  }
   else {
    mflops = 0.0;
  }
  c_print_results("FT",class,128,128,32,niter,nthreads,total_time,mflops,"          floating point",verified,"3.0 structured","25 Jun 2025","gcc","gcc","-lm","-I../common","-O3 -fopenmp -fopt-info-vec-missed=vec-miss...","-fopenmp -lm","randdp");
  if (0 == 1) 
    print_timers();
}

static void evolve(dcomplex u0[32][128][128],dcomplex u1[32][128][128],int t,int indexmap[32][128][128],int d[3])
{
  int i;
  int j;
  int k;
  
#pragma omp parallel for private (i,j,k)
  for (k = 0; k <= d[2] - 1; k += 1) {
    
#pragma omp parallel for private (i,j)
    for (j = 0; j <= d[1] - 1; j += 1) {
      
#pragma omp parallel for private (i) firstprivate (t)
      for (i = 0; i <= d[0] - 1; i += 1) {
        (u1[k][j][i] . real = u0[k][j][i] . real * ex[t * indexmap[k][j][i]] , u1[k][j][i] . imag = u0[k][j][i] . imag * ex[t * indexmap[k][j][i]]);
      }
    }
  }
}

static void compute_initial_conditions(dcomplex u0[32][128][128],int d[3])
{
  int k;
  double x0;
  double start;
  double an;
  double dummy;
  static double tmp[32769];
  int i;
  int j;
  int t;
  start = 314159265.0;
  ipow46(1220703125.0,(zstart[0] - 1) * 2 * 128 * 128 + (ystart[0] - 1) * 2 * 128,&an);
  dummy = randlc(&start,an);
  ipow46(1220703125.0,2 * 128 * 128,&an);
  for (k = 0; k <= dims[0][2] - 1; k += 1) {
    x0 = start;
    vranlc(2 * 128 * dims[0][1],&x0,1220703125.0,tmp);
    t = 1;
    
#pragma omp parallel for private (t,i,j)
    for (j = 0; j <= dims[0][1] - 1; j += 1) {
      for (i = 0; i <= 127; i += 1) {
        u0[k][j][i] . real = tmp[t++];
        u0[k][j][i] . imag = tmp[t++];
      }
    }
    if (k != dims[0][2]) 
      dummy = randlc(&start,an);
  }
}

static void ipow46(double a,int exponent,double *result)
{
  double dummy;
  double q;
  double r;
  int n;
  int n2;
   *result = 1;
  if (exponent == 0) 
    return ;
  q = a;
  r = 1;
  n = exponent;
  while(n > 1){
    n2 = n / 2;
    if (n2 * 2 == n) {
      dummy = randlc(&q,q);
      n = n2;
    }
     else {
      dummy = randlc(&r,q);
      n = n - 1;
    }
  }
  dummy = randlc(&r,q);
   *result = r;
}

static void setup()
{
  int ierr;
  int i;
  int j;
  int fstatus;
  printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version - FT Benchmark\n\n");
  niter = 6;
  printf(" Size                : %3dx%3dx%3d\n",128,128,32);
  printf(" Iterations          :     %7d\n",niter);
  
#pragma omp parallel for private (i)
  for (i = 0; i <= 2; i += 1) {
    dims[i][0] = 128;
    dims[i][1] = 128;
    dims[i][2] = 32;
  }
  
#pragma omp parallel for private (i)
  for (i = 0; i <= 2; i += 1) {
    xstart[i] = 1;
    xend[i] = 128;
    ystart[i] = 1;
    yend[i] = 128;
    zstart[i] = 1;
    zend[i] = 32;
  }
  fftblock = 16;
  fftblockpad = 18;
  if (fftblock != 16) 
    fftblockpad = fftblock + 3;
}

static void compute_indexmap(int indexmap[32][128][128],int d[3])
{
  int i;
  int j;
  int k;
  int ii;
  int ii2;
  int jj;
  int ij2;
  int kk;
  double ap;
  
#pragma omp parallel for private (ii,ii2,jj,ij2,kk,i,j,k)
  for (i = 0; i <= dims[2][0] - 1; i += 1) {
    ii = (i + 1 + xstart[2] - 2 + 128 / 2) % 128 - 128 / 2;
    ii2 = ii * ii;
    
#pragma omp parallel for private (jj,ij2,kk,j,k) firstprivate (ii2)
    for (j = 0; j <= dims[2][1] - 1; j += 1) {
      jj = (j + 1 + ystart[2] - 2 + 128 / 2) % 128 - 128 / 2;
      ij2 = jj * jj + ii2;
      
#pragma omp parallel for private (kk,k) firstprivate (ij2)
      for (k = 0; k <= dims[2][2] - 1; k += 1) {
        kk = (k + 1 + zstart[2] - 2 + 32 / 2) % 32 - 32 / 2;
        indexmap[k][j][i] = kk * kk + ij2;
      }
    }
  }
  ap = - 4.0 * 1.0e-6 * 3.141592653589793238 * 3.141592653589793238;
  ex[0] = 1.0;
  ex[1] = exp(ap);
  for (i = 2; i <= 50688; i += 1) {
    ex[i] = ex[i - 1] * ex[1];
  }
}

static void print_timers()
{
  int i;
  char *tstrings[] = {"          total ", "          setup ", "            fft ", "         evolve ", "       checksum ", "         fftlow ", "        fftcopy "};
  for (i = 0; i <= 6; i += 1) {
    if (timer_read(i) != 0.0) {
      printf("timer %2d(%16s( :%10.6f\n",i,tstrings[i],(timer_read(i)));
    }
  }
}

static void fft(int dir,dcomplex x1[32][128][128],dcomplex x2[32][128][128])
{
  dcomplex y0[128][18];
  dcomplex y1[128][18];
  if (dir == 1) {
    cffts1(1,dims[0],x1,x1,y0,y1);
    cffts2(1,dims[1],x1,x1,y0,y1);
    cffts3(1,dims[2],x1,x2,y0,y1);
  }
   else {
    cffts3(- 1,dims[2],x1,x1,y0,y1);
    cffts2(- 1,dims[1],x1,x1,y0,y1);
    cffts1(- 1,dims[0],x1,x2,y0,y1);
  }
}

static void cffts1(int is,int d[3],dcomplex x[32][128][128],dcomplex xout[32][128][128],dcomplex y0[128][18],dcomplex y1[128][18])
{
  int logd[3];
  int i;
  int j;
  int k;
  int jj;
  for (i = 0; i <= 2; i += 1) {
    logd[i] = ilog2(d[i]);
  }
{
    dcomplex y0[128][18];
    dcomplex y1[128][18];
    for (k = 0; k <= d[2] - 1; k += 1) {
      for (jj = 0; jj <= d[1] - fftblock; jj += fftblock) {
        
#pragma omp parallel for private (i,j)
        for (j = 0; j <= fftblock - 1; j += 1) {
          
#pragma omp parallel for private (i)
          for (i = 0; i <= d[0] - 1; i += 1) {
            y0[i][j] . real = x[k][j + jj][i] . real;
            y0[i][j] . imag = x[k][j + jj][i] . imag;
          }
        }
        cfftz(is,logd[0],d[0],y0,y1);
        
#pragma omp parallel for private (i,j)
        for (j = 0; j <= fftblock - 1; j += 1) {
          
#pragma omp parallel for private (i)
          for (i = 0; i <= d[0] - 1; i += 1) {
            xout[k][j + jj][i] . real = y0[i][j] . real;
            xout[k][j + jj][i] . imag = y0[i][j] . imag;
          }
        }
      }
    }
  }
}

static void cffts2(int is,int d[3],dcomplex x[32][128][128],dcomplex xout[32][128][128],dcomplex y0[128][18],dcomplex y1[128][18])
{
  int logd[3];
  int i;
  int j;
  int k;
  int ii;
  for (i = 0; i <= 2; i += 1) {
    logd[i] = ilog2(d[i]);
  }
{
    dcomplex y0[128][18];
    dcomplex y1[128][18];
    for (k = 0; k <= d[2] - 1; k += 1) {
      for (ii = 0; ii <= d[0] - fftblock; ii += fftblock) {
        
#pragma omp parallel for private (i,j)
        for (j = 0; j <= d[1] - 1; j += 1) {
          
#pragma omp parallel for private (i)
          for (i = 0; i <= fftblock - 1; i += 1) {
            y0[j][i] . real = x[k][j][i + ii] . real;
            y0[j][i] . imag = x[k][j][i + ii] . imag;
          }
        }
        cfftz(is,logd[1],d[1],y0,y1);
        
#pragma omp parallel for private (i,j)
        for (j = 0; j <= d[1] - 1; j += 1) {
          
#pragma omp parallel for private (i)
          for (i = 0; i <= fftblock - 1; i += 1) {
            xout[k][j][i + ii] . real = y0[j][i] . real;
            xout[k][j][i + ii] . imag = y0[j][i] . imag;
          }
        }
      }
    }
  }
}

static void cffts3(int is,int d[3],dcomplex x[32][128][128],dcomplex xout[32][128][128],dcomplex y0[128][18],dcomplex y1[128][18])
{
  int logd[3];
  int i;
  int j;
  int k;
  int ii;
  for (i = 0; i <= 2; i += 1) {
    logd[i] = ilog2(d[i]);
  }
{
    dcomplex y0[128][18];
    dcomplex y1[128][18];
    for (j = 0; j <= d[1] - 1; j += 1) {
      for (ii = 0; ii <= d[0] - fftblock; ii += fftblock) {
        
#pragma omp parallel for private (i,k)
        for (k = 0; k <= d[2] - 1; k += 1) {
          
#pragma omp parallel for private (i)
          for (i = 0; i <= fftblock - 1; i += 1) {
            y0[k][i] . real = x[k][j][i + ii] . real;
            y0[k][i] . imag = x[k][j][i + ii] . imag;
          }
        }
        cfftz(is,logd[2],d[2],y0,y1);
        
#pragma omp parallel for private (i,k)
        for (k = 0; k <= d[2] - 1; k += 1) {
          
#pragma omp parallel for private (i)
          for (i = 0; i <= fftblock - 1; i += 1) {
            xout[k][j][i + ii] . real = y0[k][i] . real;
            xout[k][j][i + ii] . imag = y0[k][i] . imag;
          }
        }
      }
    }
  }
}

static void fft_init(int n)
{
  int m;
  int nu;
  int ku;
  int i;
  int j;
  int ln;
  double t;
  double ti;
  nu = n;
  m = ilog2(n);
  u[0] . real = ((double )m);
  u[0] . imag = 0.0;
  ku = 1;
  ln = 1;
  for (j = 1; j <= m; j += 1) {
    t = 3.141592653589793238 / ln;
    for (i = 0; i <= ln - 1; i += 1) {
      ti = i * t;
      u[i + ku] . real = cos(ti);
      u[i + ku] . imag = sin(ti);
    }
    ku = ku + ln;
    ln = 2 * ln;
  }
}

static void cfftz(int is,int m,int n,dcomplex x[128][18],dcomplex y[128][18])
{
  int i;
  int j;
  int l;
  int mx;
  mx = ((int )u[0] . real);
  if (is != 1 && is != - 1 || m < 1 || m > mx) {
    printf("CFFTZ: Either U has not been initialized, or else\none of the input parameters is invalid%5d%5d%5d\n",is,m,mx);
    exit(1);
  }
  for (l = 1; l <= m; l += 2) {
    fftz2(is,l,m,n,fftblock,fftblockpad,u,x,y);
    if (l == m) 
      break; 
    fftz2(is,l + 1,m,n,fftblock,fftblockpad,u,y,x);
  }
  if (m % 2 == 1) {
    
#pragma omp parallel for private (i,j) firstprivate (fftblock,n)
    for (j = 0; j <= n - 1; j += 1) {
      
#pragma omp parallel for private (i)
      for (i = 0; i <= fftblock - 1; i += 1) {
        x[j][i] . real = y[j][i] . real;
        x[j][i] . imag = y[j][i] . imag;
      }
    }
  }
}

static void fftz2(int is,int l,int m,int n,int ny,int ny1,dcomplex u[128],dcomplex x[128][18],dcomplex y[128][18])
{
  int k;
  int n1;
  int li;
  int lj;
  int lk;
  int ku;
  int i;
  int j;
  int i11;
  int i12;
  int i21;
  int i22;
  dcomplex u1;
  dcomplex x11;
  dcomplex x21;
  n1 = n / 2;
  if (l - 1 == 0) {
    lk = 1;
  }
   else {
    lk = 2 << l - 1 - 1;
  }
  if (m - l == 0) {
    li = 1;
  }
   else {
    li = 2 << m - l - 1;
  }
  lj = 2 * lk;
  ku = li;
  for (i = 0; i <= li - 1; i += 1) {
    i11 = i * lk;
    i12 = i11 + n1;
    i21 = i * lj;
    i22 = i21 + lk;
    if (is >= 1) {
      u1 . real = u[ku + i] . real;
      u1 . imag = u[ku + i] . imag;
    }
     else {
      u1 . real = u[ku + i] . real;
      u1 . imag = -u[ku + i] . imag;
    }
    for (k = 0; k <= lk - 1; k += 1) {
      
#pragma omp parallel for private (j) firstprivate (i11,i12,i21,i22)
      for (j = 0; j <= ny - 1; j += 1) {
        double x11real;
        double x11imag;
        double x21real;
        double x21imag;
        x11real = x[i11 + k][j] . real;
        x11imag = x[i11 + k][j] . imag;
        x21real = x[i12 + k][j] . real;
        x21imag = x[i12 + k][j] . imag;
        y[i21 + k][j] . real = x11real + x21real;
        y[i21 + k][j] . imag = x11imag + x21imag;
        y[i22 + k][j] . real = u1 . real * (x11real - x21real) - u1 . imag * (x11imag - x21imag);
        y[i22 + k][j] . imag = u1 . real * (x11imag - x21imag) + u1 . imag * (x11real - x21real);
      }
    }
  }
}

static int ilog2(int n)
{
  int nn;
  int lg;
  if (n == 1) {
    return 0;
  }
  lg = 1;
  nn = 2;
  while(nn < n){
    nn = nn << 1;
    lg++;
  }
  return lg;
}

static void checksum(int i,dcomplex u1[32][128][128],int d[3])
{
{
    int j;
    int q;
    int r;
    int s;
    int ierr;
    dcomplex chk;
    dcomplex allchk;
    chk . real = 0.0;
    chk . imag = 0.0;
    for (j = 1; j <= 1024; j += 1) {
      q = j % 128 + 1;
      if (q >= xstart[0] && q <= xend[0]) {
        r = 3 * j % 128 + 1;
        if (r >= ystart[0] && r <= yend[0]) {
          s = 5 * j % 32 + 1;
          if (s >= zstart[0] && s <= zend[0]) {
            (chk . real = chk . real + u1[s - zstart[0]][r - ystart[0]][q - xstart[0]] . real , chk . imag = chk . imag + u1[s - zstart[0]][r - ystart[0]][q - xstart[0]] . imag);
          }
        }
      }
    }
{
      sums[i] . real += chk . real;
      sums[i] . imag += chk . imag;
    }
{
      sums[i] . real = sums[i] . real / ((double )524288);
      sums[i] . imag = sums[i] . imag / ((double )524288);
      printf("T = %5d     Checksum = %22.12e %22.12e\n",i,sums[i] . real,sums[i] . imag);
    }
  }
}

static void verify(int d1,int d2,int d3,int nt,boolean *verified,char *class)
{
  int ierr;
  int size;
  int i;
  double err;
  double epsilon;
  double vdata_real_s[7] = {(0.0), (5.546087004964e+02), (5.546385409189e+02), (5.546148406171e+02), (5.545423607415e+02), (5.544255039624e+02), (5.542683411902e+02)};
  double vdata_imag_s[7] = {(0.0), (4.845363331978e+02), (4.865304269511e+02), (4.883910722336e+02), (4.901273169046e+02), (4.917475857993e+02), (4.932597244941e+02)};
  double vdata_real_w[7] = {(0.0), (5.673612178944e+02), (5.631436885271e+02), (5.594024089970e+02), (5.560698047020e+02), (5.530898991250e+02), (5.504159734538e+02)};
  double vdata_imag_w[7] = {(0.0), (5.293246849175e+02), (5.282149986629e+02), (5.270996558037e+02), (5.260027904925e+02), (5.249400845633e+02), (5.239212247086e+02)};
  double vdata_real_a[7] = {(0.0), (5.046735008193e+02), (5.059412319734e+02), (5.069376896287e+02), (5.077892868474e+02), (5.085233095391e+02), (5.091487099959e+02)};
  double vdata_imag_a[7] = {(0.0), (5.114047905510e+02), (5.098809666433e+02), (5.098144042213e+02), (5.101336130759e+02), (5.104914655194e+02), (5.107917842803e+02)};
  double vdata_real_b[21] = {(0.0), (5.177643571579e+02), (5.154521291263e+02), (5.146409228649e+02), (5.142378756213e+02), (5.139626667737e+02), (5.137423460082e+02), (5.135547056878e+02), (5.133910925466e+02), (5.132470705390e+02), (5.131197729984e+02), (5.130070319283e+02), (5.129070537032e+02), (5.128182883502e+02), (5.127393733383e+02), (5.126691062020e+02), (5.126064276004e+02), (5.125504076570e+02), (5.125002331720e+02), (5.124551951846e+02), (5.124146770029e+02)};
  double vdata_imag_b[21] = {(0.0), (5.077803458597e+02), (5.088249431599e+02), (5.096208912659e+02), (5.101023387619e+02), (5.103976610617e+02), (5.105948019802e+02), (5.107404165783e+02), (5.108576573661e+02), (5.109577278523e+02), (5.110460304483e+02), (5.111252433800e+02), (5.111968077718e+02), (5.112616233064e+02), (5.113203605551e+02), (5.113735928093e+02), (5.114218460548e+02), (5.114656139760e+02), (5.115053595966e+02), (5.115415130407e+02), (5.115744692211e+02)};
  double vdata_real_c[21] = {(0.0), (5.195078707457e+02), (5.155422171134e+02), (5.144678022222e+02), (5.140150594328e+02), (5.137550426810e+02), (5.135811056728e+02), (5.134569343165e+02), (5.133651975661e+02), (5.132955192805e+02), (5.132410471738e+02), (5.131971141679e+02), (5.131605205716e+02), (5.131290734194e+02), (5.131012720314e+02), (5.130760908195e+02), (5.130528295923e+02), (5.130310107773e+02), (5.130103090133e+02), (5.129905029333e+02), (5.129714421109e+02)};
  double vdata_imag_c[21] = {(0.0), (5.149019699238e+02), (5.127578201997e+02), (5.122251847514e+02), (5.121090289018e+02), (5.121143685824e+02), (5.121496764568e+02), (5.121870921893e+02), (5.122193250322e+02), (5.122454735794e+02), (5.122663649603e+02), (5.122830879827e+02), (5.122965869718e+02), (5.123075927445e+02), (5.123166486553e+02), (5.123241541685e+02), (5.123304037599e+02), (5.123356167976e+02), (5.123399592211e+02), (5.123435588985e+02), (5.123465164008e+02)};
  epsilon = 1.0e-12;
   *verified = 1;
   *class = 'U';
  if (d1 == 64 && d2 == 64 && d3 == 64 && nt == 6) {
     *class = 'S';
    for (i = 1; i <= nt; i += 1) {
      err = (sums[i] . real - vdata_real_s[i]) / vdata_real_s[i];
      if (fabs(err) > epsilon) {
         *verified = 0;
        break; 
      }
      err = (sums[i] . imag - vdata_imag_s[i]) / vdata_imag_s[i];
      if (fabs(err) > epsilon) {
         *verified = 0;
        break; 
      }
    }
  }
   else if (d1 == 128 && d2 == 128 && d3 == 32 && nt == 6) {
     *class = 'W';
    for (i = 1; i <= nt; i += 1) {
      err = (sums[i] . real - vdata_real_w[i]) / vdata_real_w[i];
      if (fabs(err) > epsilon) {
         *verified = 0;
        break; 
      }
      err = (sums[i] . imag - vdata_imag_w[i]) / vdata_imag_w[i];
      if (fabs(err) > epsilon) {
         *verified = 0;
        break; 
      }
    }
  }
   else if (d1 == 256 && d2 == 256 && d3 == 128 && nt == 6) {
     *class = 'A';
    for (i = 1; i <= nt; i += 1) {
      err = (sums[i] . real - vdata_real_a[i]) / vdata_real_a[i];
      if (fabs(err) > epsilon) {
         *verified = 0;
        break; 
      }
      err = (sums[i] . imag - vdata_imag_a[i]) / vdata_imag_a[i];
      if (fabs(err) > epsilon) {
         *verified = 0;
        break; 
      }
    }
  }
   else if (d1 == 512 && d2 == 256 && d3 == 256 && nt == 20) {
     *class = 'B';
    for (i = 1; i <= nt; i += 1) {
      err = (sums[i] . real - vdata_real_b[i]) / vdata_real_b[i];
      if (fabs(err) > epsilon) {
         *verified = 0;
        break; 
      }
      err = (sums[i] . imag - vdata_imag_b[i]) / vdata_imag_b[i];
      if (fabs(err) > epsilon) {
         *verified = 0;
        break; 
      }
    }
  }
   else if (d1 == 512 && d2 == 512 && d3 == 512 && nt == 20) {
     *class = 'C';
    for (i = 1; i <= nt; i += 1) {
      err = (sums[i] . real - vdata_real_c[i]) / vdata_real_c[i];
      if (fabs(err) > epsilon) {
         *verified = 0;
        break; 
      }
      err = (sums[i] . imag - vdata_imag_c[i]) / vdata_imag_c[i];
      if (fabs(err) > epsilon) {
         *verified = 0;
        break; 
      }
    }
  }
  if (( *class) != 'U') {
    printf("Result verification successful\n");
  }
   else {
    printf("Result verification failed\n");
  }
  printf("class = %1c\n",( *class));
}
