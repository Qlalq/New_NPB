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

static double x[2*NK];
static double q[NQ];


void auto_loop_1(char size[]) {
        int j;
    for (j = 13; j >= 1; j--) {
	if (size[j] == '.') size[j] = ' ';
    }
}
void auto_loop_2(void) {
        int i;
    for (i = 0; i < 2*NK; i++) x[i] = -1.0e99;
}
void auto_loop_3(double *t1_ptr, double *t2_ptr) {
        int i;
    for ( i = 1; i <= MK+1; i++) {
	(*t2_ptr) = randlc(&(*t1_ptr), (*t1_ptr));
    }
}
void auto_loop_4(void) {
        int i;
    for ( i = 0; i <= NQ - 1; i++) {
	q[i] = 0.0;
    }
}
void auto_loop_5(double qq[]) {
        int i;
    for (i = 0; i < NQ; i++) qq[i] = 0.0;
}
void auto_loop_6(double an, int k_offset, int np, 
                           double qq[NQ], double *sx_out, double *sy_out) 
{
    // 修正点1：将所有私有变量在并行区域外、函数开头提前声明
    // 这样 OpenMP 的 pragma 指令就能“认识”它们了
    int k, i, ik, kk, l;
    double t1, t2, t3, t4, x1, x2;
    
    // 我们将 sx 和 sy 的累加结果先放在局部变量里，这是 reduction 的标准做法
    double sx = 0.0;
    double sy = 0.0;

    #pragma omp parallel for \
        private(k, i, ik, kk, l, t1, t2, t3, t4, x1, x2) \
        reduction(+:sx, sy) \
        reduction(+:qq[0:NQ])
    for (k = 1; k <= np; k++) {
        
        // 修正点2 (更安全的做法): 使用动态内存分配，避免栈溢出
        // 不再使用 double x_private[2*NK]; 这种有风险的方式
        // 而是从“堆”（Heap）上申请一块内存，这里空间大得多
        double *x_private = (double *)malloc(sizeof(double) * 2 * NK);

        // --- 接下来的计算逻辑和之前完全一样 ---

        kk = k_offset + k;
        t1 = S;
        t2 = an;

        for (i = 1; i <= 100; i++) {
            ik = kk / 2;
            if (2 * ik != kk) t3 = randlc(&t1, t2);
            if (ik == 0) break;
            t3 = randlc(&t2, t2);
            kk = ik;
        }
        
        vranlc(2*NK, &t1, A, x_private-1);

        for (i = 0; i < NK; i++) {
            x1 = 2.0 * x_private[2*i] - 1.0;
            x2 = 2.0 * x_private[2*i+1] - 1.0;
            t1 = x1*x1 + x2*x2;
            if (t1 <= 1.0) {
                t2 = sqrt(-2.0 * log(t1) / t1);
                t3 = (x1 * t2);                
                t4 = (x2 * t2);                
                l = max(fabs(t3), fabs(t4));
                qq[l] += 1.0;                
                sx = sx + t3;                
                sy = sy + t4;                
            }
        }
        
        // 修正点3：用完后释放内存，好习惯！
        free(x_private);

    } // 并行循环结束

    // 循环结束后，将所有线程累加的结果，统一加到最终的输出变量上
    *sx_out += sx;
    *sy_out += sy;
}
void auto_loop_7(double qq[]) {
          int i;
    for (i = 0; i <= NQ - 1; i++) q[i] += qq[i];
}
void auto_loop_8(double *gc_ptr) {
        int i;
    for (i = 0; i <= NQ-1; i++) {
        (*gc_ptr) = (*gc_ptr) + q[i];
    }
}

int main(int argc, char **argv) {

    double Mops, t1, t2, t3, t4, x1, x2, sx, sy, tm, an, tt, gc;
    double dum[3] = { 1.0, 1.0, 1.0 };
    int np, ierr, node, no_nodes, i, ik, kk, l, k, nit, ierrcode,
	no_large_nodes, np_add, k_offset, j;
    int nthreads = 1;
    boolean verified;
    char size[13+1];	

    sprintf(size, "%12.0f", pow(2.0, M+1));
    auto_loop_1(size);
    verified = FALSE;


    np = NN;


    vranlc(0, &(dum[0]), dum[1], &(dum[2]));
    dum[0] = randlc(&(dum[1]), dum[2]);
    
    auto_loop_2();
    
    Mops = log(sqrt(fabs(max(1.0, 1.0))));

    timer_clear(1);
    timer_clear(2);
    timer_clear(3);
    timer_start(1);

    vranlc(0, &t1, A, x);



    t1 = A;

    auto_loop_3(&t1, &t2);

    an = t1;
    tt = S;
    gc = 0.0;
    sx = 0.0;
    sy = 0.0;

    auto_loop_4();
      

    k_offset = -1;

{
    double t1, t2, t3, t4, x1, x2;
    int kk, i, ik, l;
    double qq[NQ];		

    auto_loop_5(qq);

    auto_loop_6(an, k_offset, np, qq, &sx, &sy);
    {
    auto_loop_7(qq);
    }
#if defined(_OPENMP)
    nthreads = omp_get_num_threads();
#endif     
}     

    auto_loop_8(&gc);

    timer_stop(1);
    tm = timer_read(1);

    nit = 0;
    if (M == 24) {
	if((fabs((sx- (-3.247834652034740e3))/sx) <= EPSILON) &&
	   (fabs((sy- (-6.958407078382297e3))/sy) <= EPSILON)) {
	    verified = TRUE;
	}
    } else if (M == 25) {
	if ((fabs((sx- (-2.863319731645753e3))/sx) <= EPSILON) &&
	    (fabs((sy- (-6.320053679109499e3))/sy) <= EPSILON)) {
	    verified = TRUE;
	}
    } else if (M == 28) {
	if ((fabs((sx- (-4.295875165629892e3))/sx) <= EPSILON) &&
	    (fabs((sy- (-1.580732573678431e4))/sy) <= EPSILON)) {
	    verified = TRUE;
	}
    } else if (M == 30) {
	if ((fabs((sx- (4.033815542441498e4))/sx) <= EPSILON) &&
	    (fabs((sy- (-2.660669192809235e4))/sy) <= EPSILON)) {
	    verified = TRUE;
	}
    } else if (M == 32) {
	if ((fabs((sx- (4.764367927995374e4))/sx) <= EPSILON) &&
	    (fabs((sy- (-8.084072988043731e4))/sy) <= EPSILON)) {
	    verified = TRUE;
	}
    }

    Mops = pow(2.0, M+1)/tm/1000000.0;
	  
    c_print_results("EP", CLASS, M+1, 0, 0, nit, nthreads,
		  tm, Mops, 	
		  "Random numbers generated",
		  verified, NPBVERSION, COMPILETIME,
		  CS1, CS2, CS3, CS4, CS5, CS6, CS7);

    if (TIMERS_ENABLED == TRUE) {
    }
}