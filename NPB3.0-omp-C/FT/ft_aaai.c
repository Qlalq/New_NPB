#include "npb-C.h"
#include "global.h"
static void evolve(dcomplex u0[NZ][NY][NX], dcomplex u1[NZ][NY][NX],
		   int t, int indexmap[NZ][NY][NX], int d[3]);
static void compute_initial_conditions(dcomplex u0[NZ][NY][NX], int d[3]);
static void ipow46(double a, int exponent, double *result);
static void setup(void);
static void compute_indexmap(int indexmap[NZ][NY][NX], int d[3]);
static void print_timers(void);
static void fft(int dir, dcomplex x1[NZ][NY][NX], dcomplex x2[NZ][NY][NX]);
static void cffts1(int is, int d[3], dcomplex x[NZ][NY][NX],
		   dcomplex xout[NZ][NY][NX],
		   dcomplex y0[NX][FFTBLOCKPAD],
		   dcomplex y1[NX][FFTBLOCKPAD]);
static void cffts2(int is, int d[3], dcomplex x[NZ][NY][NX],
		   dcomplex xout[NZ][NY][NX],
		   dcomplex y0[NX][FFTBLOCKPAD],
		   dcomplex y1[NX][FFTBLOCKPAD]);
static void cffts3(int is, int d[3], dcomplex x[NZ][NY][NX],
		   dcomplex xout[NZ][NY][NX],
		   dcomplex y0[NX][FFTBLOCKPAD],
		   dcomplex y1[NX][FFTBLOCKPAD]);
static void fft_init (int n);
static void cfftz (int is, int m, int n, dcomplex x[NX][FFTBLOCKPAD],
		   dcomplex y[NX][FFTBLOCKPAD]);
static void fftz2 (int is, int l, int m, int n, int ny, int ny1,
		   dcomplex u[NX], dcomplex x[NX][FFTBLOCKPAD],
		   dcomplex y[NX][FFTBLOCKPAD]);
static int ilog2(int n);
static void checksum(int i, dcomplex u1[NZ][NY][NX], int d[3]);
static void verify (int d1, int d2, int d3, int nt,
		    boolean *verified, char *class);
int main(int argc, char **argv) {
    int i, ierr;
    static dcomplex u0[NZ][NY][NX];
    static dcomplex pad1[3];
    static dcomplex u1[NZ][NY][NX];
    static dcomplex pad2[3];
    static dcomplex u2[NZ][NY][NX];
    static dcomplex pad3[3];
    static int indexmap[NZ][NY][NX];
    int iter;
    int nthreads = 1;
    double total_time, mflops;
    boolean verified;
    char class;
    for (i = 0; i < T_MAX; i++) {
	timer_clear(i);
    }
    setup();
    compute_indexmap(indexmap, dims[2]);
    compute_initial_conditions(u1, dims[0]);
    fft_init (dims[0][0]);
    fft(1, u1, u0);
    for (i = 0; i < T_MAX; i++) {
	timer_clear(i);
    }
    timer_start(T_TOTAL);
    if (TIMERS_ENABLED == TRUE) timer_start(T_SETUP);
    compute_indexmap(indexmap, dims[2]);
    compute_initial_conditions(u1, dims[0]);    
    fft_init (dims[0][0]);
    if (TIMERS_ENABLED == TRUE) {
      timer_stop(T_SETUP);
    }
    if (TIMERS_ENABLED == TRUE) {
      timer_start(T_FFT);
    }
    fft(1, u1, u0);
    if (TIMERS_ENABLED == TRUE) {   
      timer_stop(T_FFT);
    }
    for (iter = 1; iter <= niter; iter++) {
	    if (TIMERS_ENABLED == TRUE) {     
	      timer_start(T_EVOLVE);
	    }
	    evolve(u0, u1, iter, indexmap, dims[0]);
            if (TIMERS_ENABLED == TRUE) {    
	      timer_stop(T_EVOLVE);
	    }
            if (TIMERS_ENABLED == TRUE) {  
	      timer_start(T_FFT);
	    }
            fft(-1, u1, u2);
            if (TIMERS_ENABLED == TRUE) {      
	      timer_stop(T_FFT);
	    }
            if (TIMERS_ENABLED == TRUE) {     
	      timer_start(T_CHECKSUM);
	    }	
            checksum(iter, u2, dims[0]);
            if (TIMERS_ENABLED == TRUE) {   
	      timer_stop(T_CHECKSUM);
	    }
    }
    verify(NX, NY, NZ, niter, &verified, &class);
#pragma omp parallel 
  {
#if defined(_OPENMP)
#pragma omp master    
    nthreads = omp_get_num_threads();
#endif     
  } 
    timer_stop(T_TOTAL);
    total_time = timer_read(T_TOTAL);
    if( total_time != 0.0) {
	mflops = 1.0e-6*(double)(NTOTAL) *
	    (14.8157+7.19641*log((double)(NTOTAL))
	     +  (5.23518+7.21113*log((double)(NTOTAL)))*niter)
	    /total_time;
    } else {
	mflops = 0.0;
    }
    c_print_results("FT", class, NX, NY, NZ, niter, nthreads,
		    total_time, mflops, "          floating point", verified, 
		    NPBVERSION, COMPILETIME,
		    CS1, CS2, CS3, CS4, CS5, CS6, CS7);
    if (TIMERS_ENABLED == TRUE) print_timers();
}
static void evolve(dcomplex u0[NZ][NY][NX], dcomplex u1[NZ][NY][NX],
		   int t, int indexmap[NZ][NY][NX], int d[3]) {
    int i, j, k;
    for (k = 0; k < d[2]; k++) {
	for (j = 0; j < d[1]; j++) {
            for (i = 0; i < d[0]; i++) {
	      crmul(u1[k][j][i], u0[k][j][i], ex[t*indexmap[k][j][i]]);
	    }
	}
    }
}
static void compute_initial_conditions(dcomplex u0[NZ][NY][NX], int d[3]) {
    int k;
    double start, an, dummy;
    
    // 初始的种子计算保持不变
    start = SEED;
    // 下面这两行是为了计算出整个任务的“总起点”
    ipow46(A, (zstart[0]-1)*2*NX*NY + (ystart[0]-1)*2*NX, &an);
    dummy = randlc(&start, an);
    
    // 就像召集了很多同学，准备分工合作
    // schedule(static) 表示平均分配任务，每个人拿到的任务量差不多
    #pragma omp parallel for private(k, an, dummy) schedule(static)
    for (k = 0; k < dims[0][2]; k++) {
        // --- 这是每个同学都需要独立拥有的“工具”和“草稿纸” ---
        double x0;
        // tmp 数组很大，放在每个线程自己的栈上，避免互相干扰
        double tmp[NX*2*MAXDIM+1]; 
        int i, j, t;

        // --- 每个同学计算自己第一个任务需要的“接力棒” ---
        // ipow46 和 randlc 不是线程安全的，所以这里的 an 和 dummy 必须是私有的
        // 我们根据 k 的值，从“总起点”start 出发，计算出当前 k 需要的精确起点
        x0 = start;
        ipow46(A, (double)k * 2.0 * NX * dims[0][1], &an);
        dummy = randlc(&x0, an);

        // --- 开始独立工作！---
        // 用自己算好的 x0 作为种子，生成随机数
        vranlc(2*NX*dims[0][1], &x0, A, tmp);

        t = 1;
        for (j = 0; j < dims[0][1]; j++) {
            for (i = 0; i < NX; i++) {
                u0[k][j][i].real = tmp[t++];
                u0[k][j][i].imag = tmp[t++];
            }
        }
        // 注意：原代码结尾的 randlc 更新全局 start 的步骤被去掉了
        // 因为我们已经打破了这种依赖关系！
    }
}
static void ipow46(double a, int exponent, double *result) {
    double dummy, q, r;
    int n, n2;
    *result = 1;
    if (exponent == 0) return;
    q = a;
    r = 1;
    n = exponent;
    while (n > 1) {
	n2 = n/2;
	if (n2 * 2 == n) {
            dummy = randlc(&q, q);
            n = n2;
	} else {
            dummy = randlc(&r, q);
            n = n-1;
	}
    }
    dummy = randlc(&r, q);
    *result = r;
}
static void setup(void) {
    int ierr, i, j, fstatus;
    printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version"
	   " - FT Benchmark\n\n");
    niter = NITER_DEFAULT;
    printf(" Size                : %3dx%3dx%3d\n", NX, NY, NZ);
    printf(" Iterations          :     %7d\n", niter);
    for (i = 0;i < 3 ; i++) {
	dims[i][0] = NX;
	dims[i][1] = NY;
	dims[i][2] = NZ;
    }
    for (i = 0; i < 3; i++) {
	xstart[i] = 1;
	xend[i]   = NX;
	ystart[i] = 1;
        yend[i]   = NY;
        zstart[i] = 1;
        zend[i]   = NZ;
    }
    fftblock = FFTBLOCK_DEFAULT;
    fftblockpad = FFTBLOCKPAD_DEFAULT;
    if (fftblock != FFTBLOCK_DEFAULT) fftblockpad = fftblock+3;
}
static void compute_indexmap(int indexmap[NZ][NY][NX], int d[3]) {
    int i, j, k, ii, ii2, jj, ij2, kk;
    double ap;
    for (i = 0; i < dims[2][0]; i++) {
	ii =  (i+1+xstart[2]-2+NX/2)%NX - NX/2;
	ii2 = ii*ii;
	for (j = 0; j < dims[2][1]; j++) {
            jj = (j+1+ystart[2]-2+NY/2)%NY - NY/2;
            ij2 = jj*jj+ii2;
            for (k = 0; k < dims[2][2]; k++) {
		kk = (k+1+zstart[2]-2+NZ/2)%NZ - NZ/2;
		indexmap[k][j][i] = kk*kk+ij2;
	    }
	}
    }
    ap = - 4.0 * ALPHA * PI * PI;
    ex[0] = 1.0;
    ex[1] = exp(ap);
    for (i = 2; i <= EXPMAX; i++) {
	ex[i] = ex[i-1]*ex[1];
    }
}
static void print_timers(void) {
    int i;
    char *tstrings[] = { "          total ",
			 "          setup ", 
			 "            fft ", 
			 "         evolve ", 
			 "       checksum ", 
			 "         fftlow ", 
			 "        fftcopy " };
    for (i = 0; i < T_MAX; i++) {
	if (timer_read(i) != 0.0) {
            printf("timer %2d(%16s( :%10.6f\n", i, tstrings[i], timer_read(i));
	}
    }
}
static void fft(int dir, dcomplex x1[NZ][NY][NX], dcomplex x2[NZ][NY][NX]) {
    dcomplex y0[NX][FFTBLOCKPAD];
    dcomplex y1[NX][FFTBLOCKPAD];
    if (dir == 1) {
        cffts1(1, dims[0], x1, x1, y0, y1);	
        cffts2(1, dims[1], x1, x1, y0, y1);	
        cffts3(1, dims[2], x1, x2, y0, y1);	
    } else {
	cffts3(-1, dims[2], x1, x1, y0, y1);	
    cffts2(-1, dims[1], x1, x1, y0, y1);	
    cffts1(-1, dims[0], x1, x2, y0, y1);	
    }
}
static void cffts1(int is, int d[3], dcomplex x[NZ][NY][NX],
		   dcomplex xout[NZ][NY][NX],
		   dcomplex y0_dummy[NX][FFTBLOCKPAD], // 这两个参数我们依然保留，以防其他地方调用时出错
		   dcomplex y1_dummy[NX][FFTBLOCKPAD]) {
    int logd[3];
    int i, j, k, jj;

    for (i = 0; i < 3; i++) {
        logd[i] = ilog2(d[i]);
    }

    // 魔法咒语改动在这里！我们把 y0 和 y1 从 private 子句中移除了
    // 因为它们在循环内部定义，所以已经是私有的了
    #pragma omp parallel for private(k, jj, j, i) schedule(static)
    for (k = 0; k < d[2]; k++) {
        // --- 这两个“工作台”在并行循环内定义，每个线程都会有自己的一份 ---
        dcomplex y0[NX][FFTBLOCKPAD];
        dcomplex y1[NX][FFTBLOCKPAD];
        
        for (jj = 0; jj <= d[1] - fftblock; jj+=fftblock) {
            // 1. 把数据从大仓库 x 搬到自己的工作台 y0
            for (j = 0; j < fftblock; j++) {
                for (i = 0; i < d[0]; i++) {
                    y0[i][j].real = x[k][j+jj][i].real;
                    y0[i][j].imag = x[k][j+jj][i].imag;
                }
            }

            // 2. 调用“高级技工” cfftz 加工自己工作台上的数据
            cfftz (is, logd[0], d[0], y0, y1);

            // 3. 把加工好的数据从工作台 y0 搬回大仓库 xout
            for (j = 0; j < fftblock; j++) {
                for (i = 0; i < d[0]; i++) {
                    xout[k][j+jj][i].real = y0[i][j].real;
                    xout[k][j+jj][i].imag = y0[i][j].imag;
                }
            }
        }
    }
}
static void cffts2(int is, int d[3], dcomplex x[NZ][NY][NX],
		   dcomplex xout[NZ][NY][NX],
		   dcomplex y0_dummy[NX][FFTBLOCKPAD], // 保留参数，保持函数接口一致
		   dcomplex y1_dummy[NX][FFTBLOCKPAD]) {
    int logd[3];
    int i, j, k, ii;

    for (i = 0; i < 3; i++) {
        logd[i] = ilog2(d[i]);
    }

    // 和 cffts1 完全相同的优化方法！
    // 只是里面的循环变量和数组索引有些许不同
    #pragma omp parallel for private(k, ii, j, i) schedule(static)
    for (k = 0; k < d[2]; k++) {
        // 每个线程都有自己私有的 y0 和 y1，互不干扰
        dcomplex y0[NX][FFTBLOCKPAD];
        dcomplex y1[NX][FFTBLOCKPAD];

        for (ii = 0; ii <= d[0] - fftblock; ii+=fftblock) {
            // 1. 数据搬运：从大仓库 x -> 私人工作台 y0
	    for (j = 0; j < d[1]; j++) {
		for (i = 0; i < fftblock; i++) {
		    y0[j][i].real = x[k][j][i+ii].real;
		    y0[j][i].imag = x[k][j][i+ii].imag;
		}
	    }
            
            // 2. 核心加工：调用 cfftz
	    cfftz (is, logd[1], d[1], y0, y1);

            // 3. 结果搬回：从私人工作台 y0 -> 大仓库 xout
            for (j = 0; j < d[1]; j++) {
                for (i = 0; i < fftblock; i++) {
                    xout[k][j][i+ii].real = y0[j][i].real;
                    xout[k][j][i+ii].imag = y0[j][i].imag;
		}
	    }
	}
    }
}
static void cffts3(int is, int d[3], dcomplex x[NZ][NY][NX],
		   dcomplex xout[NZ][NY][NX],
		   dcomplex y0_dummy[NX][FFTBLOCKPAD], // 保持函数接口一致
		   dcomplex y1_dummy[NX][FFTBLOCKPAD]) {
    int logd[3];
    int i, j, k, ii;

    for (i = 0; i < 3; i++) {
        logd[i] = ilog2(d[i]);
    }

    // 魔法咒语！这次我们并行化 j 循环
    #pragma omp parallel for private(j, ii, k, i) schedule(static)
    for (j = 0; j < d[1]; j++) {
        // 每个线程都在自己的空间里创建私有的工作台
        dcomplex y0[NX][FFTBLOCKPAD];
        dcomplex y1[NX][FFTBLOCKPAD];

        for (ii = 0; ii <= d[0] - fftblock; ii+=fftblock) {
            // 1. 数据搬运
	    for (k = 0; k < d[2]; k++) {
		for (i = 0; i < fftblock; i++) {
		    y0[k][i].real = x[k][j][i+ii].real;
		    y0[k][i].imag = x[k][j][i+ii].imag;
		}
	    }
            
            // 2. 核心加工
            cfftz (is, logd[2], d[2], y0, y1);
            
            // 3. 结果搬回
            for (k = 0; k < d[2]; k++) {
                for (i = 0; i < fftblock; i++) {
                    xout[k][j][i+ii].real = y0[k][i].real;
                    xout[k][j][i+ii].imag = y0[k][i].imag;
		}
	    }
	}
    }
}
static void fft_init (int n) {
    int m,nu,ku,i,j,ln;
    double t, ti;
    nu = n;
    m = ilog2(n);
    u[0].real = (double)m;
    u[0].imag = 0.0;
    ku = 1;
    ln = 1;
    for (j = 1; j <= m; j++) {
	t = PI / ln;
	for (i = 0; i <= ln - 1; i++) {
            ti = i * t;
            u[i+ku].real = cos(ti);
	    u[i+ku].imag = sin(ti);
	}
	ku = ku + ln;
	ln = 2 * ln;
    }
}
static void cfftz (int is, int m, int n, dcomplex x[NX][FFTBLOCKPAD],
		   dcomplex y[NX][FFTBLOCKPAD]) {
    int i,j,l,mx;
    mx = (int)(u[0].real);
    if ((is != 1 && is != -1) || m < 1 || m > mx) {
	printf("CFFTZ: Either U has not been initialized, or else\n"
	       "one of the input parameters is invalid%5d%5d%5d\n",
	       is, m, mx);
	exit(1);
    }
    for (l = 1; l <= m; l+=2) {
        fftz2 (is, l, m, n, fftblock, fftblockpad, u, x, y);
        if (l == m) break;
	fftz2 (is, l + 1, m, n, fftblock, fftblockpad, u, y, x);
    }
    if (m % 2 == 1) {
	for (j = 0; j < n; j++) {
	    for (i = 0; i < fftblock; i++) {
		x[j][i].real = y[j][i].real;
		x[j][i].imag = y[j][i].imag;
	    }
	}
    }
}
static void fftz2 (int is, int l, int m, int n, int ny, int ny1,
		   dcomplex u[NX], dcomplex x[NX][FFTBLOCKPAD],
		   dcomplex y[NX][FFTBLOCKPAD]) {
    int k,n1,li,lj,lk,ku,i,j,i11,i12,i21,i22;
    dcomplex u1,x11,x21;
    n1 = n / 2;
    if (l-1 == 0) {
	lk = 1;
    } else {
	lk = 2 << ((l - 1)-1);
    }
    if (m-l == 0) {
	li = 1;
    } else {
	li = 2 << ((m - l)-1);
    }
    lj = 2 * lk;
    ku = li;
    for (i = 0; i < li; i++) {
        i11 = i * lk;
        i12 = i11 + n1;
        i21 = i * lj;
        i22 = i21 + lk;
        if (is >= 1) {
          u1.real = u[ku+i].real;
          u1.imag = u[ku+i].imag;
        } else {
          u1.real = u[ku+i].real;
          u1.imag = -u[ku+i].imag;
        }
        for (k = 0; k < lk; k++) {
	    for (j = 0; j < ny; j++) {
		double x11real, x11imag;
		double x21real, x21imag;
		x11real = x[i11+k][j].real;
		x11imag = x[i11+k][j].imag;
		x21real = x[i12+k][j].real;
		x21imag = x[i12+k][j].imag;
		y[i21+k][j].real = x11real + x21real;
		y[i21+k][j].imag = x11imag + x21imag;
		y[i22+k][j].real = u1.real * (x11real - x21real)
		    - u1.imag * (x11imag - x21imag);
		y[i22+k][j].imag = u1.real * (x11imag - x21imag)
		    + u1.imag * (x11real - x21real);
	    }
	}
    }
}
static int ilog2(int n) {
    int nn, lg;
    if (n == 1) {
	return 0;
    }
    lg = 1;
    nn = 2;
    while (nn < n) {
	nn = nn << 1;
	lg++;
    }
    return lg;
}
static void checksum(int i, dcomplex u1[NZ][NY][NX], int d[3]) {
{
    int j, q,r,s, ierr;
    dcomplex chk,allchk;
    chk.real = 0.0;
    chk.imag = 0.0;
    for (j = 1; j <= 1024; j++) {
	q = j%NX+1;
	if (q >= xstart[0] && q <= xend[0]) {
            r = (3*j)%NY+1;
            if (r >= ystart[0] && r <= yend[0]) {
		s = (5*j)%NZ+1;
		if (s >= zstart[0] && s <= zend[0]) {
		  cadd(chk,chk,u1[s-zstart[0]][r-ystart[0]][q-xstart[0]]);
		}
	    }
	}
    }
    {
	sums[i].real += chk.real;
	sums[i].imag += chk.imag;
    }
  {    
    sums[i].real = sums[i].real/(double)(NTOTAL);
    sums[i].imag = sums[i].imag/(double)(NTOTAL);
    printf("T = %5d     Checksum = %22.12e %22.12e\n",
	   i, sums[i].real, sums[i].imag);
  }
}
}
static void verify (int d1, int d2, int d3, int nt,
		    boolean *verified, char *class) {
    int ierr, size, i;
    double err, epsilon;
    double vdata_real_s[6+1] = { 0.0,
				 5.546087004964e+02,
				 5.546385409189e+02,
				 5.546148406171e+02,
				 5.545423607415e+02,
				 5.544255039624e+02,
				 5.542683411902e+02 };
    double vdata_imag_s[6+1] = { 0.0,
				 4.845363331978e+02,
				 4.865304269511e+02,
				 4.883910722336e+02,
				 4.901273169046e+02,
				 4.917475857993e+02,
				 4.932597244941e+02 };
    double vdata_real_w[6+1] = { 0.0,
				 5.673612178944e+02,
				 5.631436885271e+02,
				 5.594024089970e+02,
				 5.560698047020e+02,
				 5.530898991250e+02,
				 5.504159734538e+02 };
    double vdata_imag_w[6+1] = { 0.0,
				 5.293246849175e+02,
				 5.282149986629e+02,
				 5.270996558037e+02, 
				 5.260027904925e+02, 
				 5.249400845633e+02,
				 5.239212247086e+02 };
    double vdata_real_a[6+1] = { 0.0,
				 5.046735008193e+02,
				 5.059412319734e+02,
				 5.069376896287e+02,
				 5.077892868474e+02,
				 5.085233095391e+02,
				 5.091487099959e+02 };
    double vdata_imag_a[6+1] = { 0.0,
				 5.114047905510e+02,
				 5.098809666433e+02,
				 5.098144042213e+02,
				 5.101336130759e+02,
				 5.104914655194e+02,
				 5.107917842803e+02 };
    double vdata_real_b[20+1] = { 0.0,
				  5.177643571579e+02,
				  5.154521291263e+02,
				  5.146409228649e+02,
				  5.142378756213e+02,
				  5.139626667737e+02,
				  5.137423460082e+02,
				  5.135547056878e+02,
				  5.133910925466e+02,
				  5.132470705390e+02,
				  5.131197729984e+02,
				  5.130070319283e+02,
				  5.129070537032e+02,
				  5.128182883502e+02,
				  5.127393733383e+02,
				  5.126691062020e+02,
				  5.126064276004e+02,
				  5.125504076570e+02,
				  5.125002331720e+02,
				  5.124551951846e+02,
				  5.124146770029e+02 };
    double vdata_imag_b[20+1] = { 0.0,
				  5.077803458597e+02,
				  5.088249431599e+02,                  
				  5.096208912659e+02,
				  5.101023387619e+02,                  
				  5.103976610617e+02,                  
				  5.105948019802e+02,                  
				  5.107404165783e+02,                  
				  5.108576573661e+02,                  
				  5.109577278523e+02,
				  5.110460304483e+02,                  
				  5.111252433800e+02,                  
				  5.111968077718e+02,                  
				  5.112616233064e+02,                  
				  5.113203605551e+02,                  
				  5.113735928093e+02,                  
				  5.114218460548e+02,
				  5.114656139760e+02,
				  5.115053595966e+02,
				  5.115415130407e+02,
				  5.115744692211e+02 };
    double vdata_real_c[20+1] = { 0.0,
				  5.195078707457e+02,
				  5.155422171134e+02,
				  5.144678022222e+02,
				  5.140150594328e+02,
				  5.137550426810e+02,
				  5.135811056728e+02,
				  5.134569343165e+02,
				  5.133651975661e+02,
				  5.132955192805e+02,
				  5.132410471738e+02,
				  5.131971141679e+02,
				  5.131605205716e+02,
				  5.131290734194e+02,
				  5.131012720314e+02,
				  5.130760908195e+02,
				  5.130528295923e+02,
				  5.130310107773e+02,
				  5.130103090133e+02,
				  5.129905029333e+02,
				  5.129714421109e+02 };
    double vdata_imag_c[20+1] = { 0.0,
				  5.149019699238e+02,
				  5.127578201997e+02,
				  5.122251847514e+02,
				  5.121090289018e+02,
				  5.121143685824e+02,
				  5.121496764568e+02,
				  5.121870921893e+02,
				  5.122193250322e+02,
				  5.122454735794e+02,
				  5.122663649603e+02,
				  5.122830879827e+02,
				  5.122965869718e+02,
				  5.123075927445e+02,
				  5.123166486553e+02,
				  5.123241541685e+02,
				  5.123304037599e+02,
				  5.123356167976e+02,
				  5.123399592211e+02,
				  5.123435588985e+02,
				  5.123465164008e+02 };
    epsilon = 1.0e-12;
    *verified = TRUE;
    *class = 'U';
    if (d1 == 64 &&
	d2 == 64 &&
	d3 == 64 &&
	nt == 6) {
	*class = 'S';
	for (i = 1; i <= nt; i++) {
            err = (get_real(sums[i]) - vdata_real_s[i]) / vdata_real_s[i];
            if (fabs(err) > epsilon) {
	      *verified = FALSE;
	      break;
	    }
            err = (get_imag(sums[i]) - vdata_imag_s[i]) / vdata_imag_s[i];
            if (fabs(err) > epsilon) {
	      *verified = FALSE;
	      break;
	    }
	}
    } else if (d1 == 128 &&
	       d2 == 128 &&
	       d3 == 32 &&
	       nt == 6) {
	*class = 'W';
	for (i = 1; i <= nt; i++) {
            err = (get_real(sums[i]) - vdata_real_w[i]) / vdata_real_w[i];
            if (fabs(err) > epsilon) {
	      *verified = FALSE;
	      break;
	    }
            err = (get_imag(sums[i]) - vdata_imag_w[i]) / vdata_imag_w[i];
            if (fabs(err) > epsilon) {
	      *verified = FALSE;
	      break;
	    }
	}
    } else if (d1 == 256 &&
	       d2 == 256 &&
	       d3 == 128 &&
	       nt == 6) {
	*class = 'A';
	for (i = 1; i <= nt; i++) {
            err = (get_real(sums[i]) - vdata_real_a[i]) / vdata_real_a[i];
            if (fabs(err) > epsilon) {
	      *verified = FALSE;
	      break;
	    }
            err = (get_imag(sums[i]) - vdata_imag_a[i]) / vdata_imag_a[i];
            if (fabs(err) > epsilon) {
	      *verified = FALSE;
	      break;
	    }
	}
    } else if (d1 == 512 &&
	       d2 == 256 &&
	       d3 == 256 &&
	       nt == 20) {
	*class = 'B';
	for (i = 1; i <= nt; i++) {
            err = (get_real(sums[i]) - vdata_real_b[i]) / vdata_real_b[i];
            if (fabs(err) > epsilon) {
	      *verified = FALSE;
	      break;
	    }
            err = (get_imag(sums[i]) - vdata_imag_b[i]) / vdata_imag_b[i];
            if (fabs(err) > epsilon) {
	      *verified = FALSE;
	      break;
	    }
	}
    } else if (d1 == 512 &&
	       d2 == 512 &&
	       d3 == 512 &&
	       nt == 20) {
	*class = 'C';
	for (i = 1; i <= nt; i++) {
            err = (get_real(sums[i]) - vdata_real_c[i]) / vdata_real_c[i];
            if (fabs(err) > epsilon) {
	      *verified = FALSE;
	      break;
	    }
            err = (get_imag(sums[i]) - vdata_imag_c[i]) / vdata_imag_c[i];
            if (fabs(err) > epsilon) {
	      *verified = FALSE;
	      break;
	    }
	}
    }
    if (*class != 'U') {
	printf("Result verification successful\n");
    } else {
	printf("Result verification failed\n");
    }
    printf("class = %1c\n", *class);
}