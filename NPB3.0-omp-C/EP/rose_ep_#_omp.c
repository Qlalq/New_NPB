#include "npb-C.h"
#include "npbparams.h"
#define	MK		16
#define	MM		(M - MK)
#define	NN		(1 << MM)
#define	NK		(1 << MK)
#define	NQ		10
#define EPSILON		1.0e-8
#define	A		1220703125.0
#define	S		271828183.0
#define	TIMERS_ENABLED	FALSE
#include <omp.h> 
static double x[131072];
static double q[10];

int main(int argc,char **argv)
{
  double Mops;
  double t1;
  double t2;
  double t3;
  double t4;
  double x1;
  double x2;
  double sx;
  double sy;
  double tm;
  double an;
  double tt;
  double gc;
  double dum[3] = {(1.0), (1.0), (1.0)};
  int np;
  int ierr;
  int node;
  int no_nodes;
  int i;
  int ik;
  int kk;
  int l;
  int k;
  int nit;
  int ierrcode;
  int no_large_nodes;
  int np_add;
  int k_offset;
  int j;
  int nthreads = 1;
  boolean verified;
  char size[14];
  sprintf(size,"%12.0f",(pow(2.0,(25 + 1))));
  
#pragma omp parallel for private (j)
  for (j = 13; j >= 1; j += -1) {
    if (size[j] == '.') 
      size[j] = ' ';
  }
  verified = 0;
  np = 1 << 25 - 16;
  vranlc(0,&dum[0],dum[1],&dum[2]);
  dum[0] = randlc(&dum[1],dum[2]);
  
#pragma omp parallel for private (i)
  for (i = 0; i <= 131071; i += 1) {
    x[i] = - 1.0e99;
  }
  Mops = log((sqrt((fabs((1.0 > 1.0?1.0 : 1.0))))));
  timer_clear(1);
  timer_clear(2);
  timer_clear(3);
  timer_start(1);
  vranlc(0,&t1,1220703125.0,x);
  t1 = 1220703125.0;
  for (i = 1; i <= 17; i += 1) {
    t2 = randlc(&t1,t1);
  }
  an = t1;
  tt = 271828183.0;
  gc = 0.0;
  sx = 0.0;
  sy = 0.0;
  
#pragma omp parallel for private (i)
  for (i = 0; i <= 9; i += 1) {
    q[i] = 0.0;
  }
  k_offset = - 1;
{
    double t1;
    double t2;
    double t3;
    double t4;
    double x1;
    double x2;
    int kk;
    int i;
    int ik;
    int l;
    double qq[10];
    
#pragma omp parallel for private (i)
    for (i = 0; i <= 9; i += 1) {
      qq[i] = 0.0;
    }
    for (k = 1; k <= np; k += 1) {
      kk = k_offset + k;
      t1 = 271828183.0;
      t2 = an;
      for (i = 1; i <= 100; i += 1) {
        ik = kk / 2;
        if (2 * ik != kk) 
          t3 = randlc(&t1,t2);
        if (ik == 0) 
          break; 
        t3 = randlc(&t2,t2);
        kk = ik;
      }
      if (0 == 1) 
        timer_start(3);
      vranlc(2 * (1 << 16),&t1,1220703125.0,x - 1);
      if (0 == 1) 
        timer_stop(3);
      if (0 == 1) 
        timer_start(2);
      for (i = 0; i <= 65535; i += 1) {
        x1 = 2.0 * x[2 * i] - 1.0;
        x2 = 2.0 * x[2 * i + 1] - 1.0;
        t1 = x1 * x1 + x2 * x2;
        if (t1 <= 1.0) {
          t2 = sqrt(- 2.0 * log(t1) / t1);
          t3 = x1 * t2;
          t4 = x2 * t2;
          l = ((fabs(t3) > fabs(t4)?fabs(t3) : fabs(t4)));
          qq[l] += 1.0;
          sx = sx + t3;
          sy = sy + t4;
        }
      }
      if (0 == 1) 
        timer_stop(2);
    }
{
      
#pragma omp parallel for private (i)
      for (i = 0; i <= 9; i += 1) {
        q[i] += qq[i];
      }
    }
#if defined(_OPENMP)
#endif     
  }
  
#pragma omp parallel for private (i) reduction (+:gc)
  for (i = 0; i <= 9; i += 1) {
    gc = gc + q[i];
  }
  timer_stop(1);
  tm = timer_read(1);
  nit = 0;
  if (25 == 24) {
    if (fabs((sx - - 3.247834652034740e3) / sx) <= 1.0e-8 && fabs((sy - - 6.958407078382297e3) / sy) <= 1.0e-8) {
      verified = 1;
    }
  }
   else if (25 == 25) {
    if (fabs((sx - - 2.863319731645753e3) / sx) <= 1.0e-8 && fabs((sy - - 6.320053679109499e3) / sy) <= 1.0e-8) {
      verified = 1;
    }
  }
   else if (25 == 28) {
    if (fabs((sx - - 4.295875165629892e3) / sx) <= 1.0e-8 && fabs((sy - - 1.580732573678431e4) / sy) <= 1.0e-8) {
      verified = 1;
    }
  }
   else if (25 == 30) {
    if (fabs((sx - 4.033815542441498e4) / sx) <= 1.0e-8 && fabs((sy - - 2.660669192809235e4) / sy) <= 1.0e-8) {
      verified = 1;
    }
  }
   else if (25 == 32) {
    if (fabs((sx - 4.764367927995374e4) / sx) <= 1.0e-8 && fabs((sy - - 8.084072988043731e4) / sy) <= 1.0e-8) {
      verified = 1;
    }
  }
  Mops = pow(2.0,(25 + 1)) / tm / 1000000.0;
  c_print_results("EP",'W',25 + 1,0,0,nit,nthreads,tm,Mops,"Random numbers generated",verified,"3.0 structured","25 Jun 2025","gcc","gcc","-lm","-I../common","-O3 -fopenmp -fopt-info-vec-missed=vec-miss...","-fopenmp -lm","randdp");
  if (0 == 1) {
  }
}
