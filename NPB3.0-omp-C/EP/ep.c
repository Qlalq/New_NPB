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

// Global arrays: x is written to repeatedly, q accumulates results
static double x[2*NK];
static double q[NQ];

int main(int argc, char **argv) {

    // Variables for timing, performance, verification
    double Mops, t1_glob, t2_glob, t3_glob, t4_glob, x1_glob, x2_glob, sx, sy, tm, an, tt, gc;
    double dum[3] = { 1.0, 1.0, 1.0 };
    int np, ierr, node, no_nodes, i_glob, ik_glob, kk_glob, l_glob, k_glob, nit, ierrcode,
	no_large_nodes, np_add, k_offset, j_glob; // Added _glob suffix to distinguish from potential loop indices
    int nthreads = 1; // Will be set by OpenMP
    boolean verified;
    char size[13+1];	

    // Initial setup and verification parameters
    sprintf(size, "%12.0f", pow(2.0, M+1));
    for (j_glob = 13; j_glob >= 1; j_glob--) {
	if (size[j_glob] == '.') size[j_glob] = ' ';
    }
    verified = FALSE;

    np = NN;

    // Initial random number generation calls - likely for setup/warmup
    vranlc(0, &(dum[0]), dum[1], &(dum[2]));
    dum[0] = randlc(&(dum[1]), dum[2]);
    
    // Initialize array x - small loop, not the main workload
    for (i_glob = 0; i_glob < 2*NK; i_glob++) x[i_glob] = -1.0e99;
    
    // Timer initializations
    Mops = log(sqrt(fabs(max(1.0, 1.0))));

    timer_clear(1);
    timer_clear(2);
    timer_clear(3);
    timer_start(1);

    // Another random number generation call - purpose unclear without vranlc source
    vranlc(0, &t1_glob, A, x);

    // Seed generation for the main loop
    t1_glob = A;
    for ( i_glob = 1; i_glob <= MK+1; i_glob++) {
	t2_glob = randlc(&t1_glob, t1_glob);
    }
    an = t1_glob;
    tt = S;
    gc = 0.0;
    sx = 0.0; // Reduction variable
    sy = 0.0; // Reduction variable

    // Initialize q array - small loop
    for ( i_glob = 0; i_glob <= NQ - 1; i_glob++) {
	q[i_glob] = 0.0; // Reduction array
    }
      
    k_offset = -1;

    // --- Main computational block ---
    // This block contains the core loops. Variables declared here (like qq)
    // are local to this block's scope.
{
    // Local variables for the main computation block.
    // t1, t2, t3, t4, x1, x2, kk, i, ik, l are used inside the loops.
    // qq is an array used for accumulating counts within this block.
    double t1, t2, t3, t4, x1, x2;
    int kk, i, ik, l;
    double qq[NQ];		// Array for accumulating counts within this block

    // Initialize the local accumulator array qq
    for (i = 0; i < NQ; i++) qq[i] = 0.0;

    // The main outer loop (k) iterates 'np' times.
    // Dependency Analysis for k loop:
    // - Each iteration k calculates a seed based on k_offset + k. This seed generation
    //   is independent across k iterations.
    // - The call vranlc(2*NK, &t1, A, x-1) generates 2*NK random numbers and writes
    //   them into the GLOBAL static array 'x'. Since the same address 'x-1' is
    //   passed every time, each iteration k *overwrites* the random numbers generated
    //   by the previous iteration in 'x'.
    // - The inner loop (i) processes the numbers *currently* in 'x'.
    // Conclusion: Due to the serial write/read pattern on the global 'x' array
    // imposed by the vranlc call within the k loop, the k loop itself cannot
    // be safely parallelized by simply distributing iterations without
    // major changes to how random numbers are generated or stored (e.g.,
    // generating into thread-private buffers). Since we are not allowed to
    // modify vranlc or its call signature, the k loop must remain serial
    // to avoid a race condition on the 'x' array writes by vranlc.

    for (k_glob = 1; k_glob <= np; k_glob++) { // This loop remains serial
        // Variables local to k iteration (some already declared in block scope)
	    kk = k_offset + k_glob; // k_glob used instead of k to differentiate from loop index

	    // Seed generation for this k iteration - Uses local t1, t2
        t1 = S;
	    t2 = an;
	    for (i = 1; i <= 100; i++) {
                ik = kk / 2;
                if (2 * ik != kk) t3 = randlc(&t1, t2);
                if (ik == 0) break;
                t3 = randlc(&t2, t2);
                kk = ik; // kk is modified, iteration-specific
	    }

        // Generate random numbers for this k iteration into the global x array
        // This must happen serially for each k due to the write dependency on x
	    if (TIMERS_ENABLED == TRUE) timer_start(3);
        // vranlc writes 2*NK numbers into the global static array x starting at x
        // This is the operation that prevents simple parallelization of the k loop
	    vranlc(2*NK, &t1, A, x-1); // t1 is the seed calculated above
	    if (TIMERS_ENABLED == TRUE) timer_stop(3);

	    if (TIMERS_ENABLED == TRUE) timer_start(2);

        // The inner loop (i) processes the numbers *just* generated in x for this k.
        // Dependency Analysis for i loop:
        // - Reads from x (read-only within this loop iteration).
        // - Calculations of x1, x2, t1, t2, t3, t4 are independent for each i.
        // - Updates to qq[l], sx, sy.
        // - qq[l] += 1.0 updates an element of the qq array. If the i loop is
        //   parallelized, multiple threads will update elements of the *same*
        //   qq array concurrently. This requires synchronization (atomic or critical).
        // - sx += t3 and sy += t4 are sum reductions. These can be handled
        //   efficiently with OpenMP reduction clauses.
        // Conclusion: The inner loop is suitable for parallelization provided
        // updates to qq are synchronized and sx/sy are handled as reductions.

        // This loop can be parallelized for different values of 'i'
	    for ( i = 0; i < NK; i++) { // This loop can be parallelized
            // Variables local to i iteration
            // x[2*i] and x[2*i+1] are read. Accesses are disjoint for different i.
            x1 = 2.0 * x[2*i] - 1.0;
            x2 = 2.0 * x[2*i+1] - 1.0;
            t1 = pow2(x1) + pow2(x2); // t1, t2, t3, t4 here are different from seed calculation
            if (t1 <= 1.0) {
		        t2 = sqrt(-2.0 * log(t1) / t1);
		        t3 = (x1 * t2);				// Part of sx sum
		        t4 = (x2 * t2);				// Part of sy sum
		        l = max(fabs(t3), fabs(t4)); // l is loop-local

                // Accumulate count in qq array. This is a write to shared memory
                // (the qq array local to the block, but shared by threads in parallel i loop).
                // Needs synchronization if the i loop is parallel.
		        qq[l] += 1.0; // Potential race if i loop is parallel

                // Accumulate sums for sx and sy. These are reductions.
		        sx = sx + t3; // Reduction operation
		        sy = sy + t4; // Reduction operation
            }
	    } // End of inner loop 'i' - parallelization potential here

	    if (TIMERS_ENABLED == TRUE) timer_stop(2);
    } // End of outer loop 'k' - serial execution due to vranlc/x dependency

    // After the k loop finishes (serially), accumulate the counts from the
    // block-local qq array into the global q array.
    // This is a simple serial loop after the main computation.
    {
      for (i = 0; i <= NQ - 1; i++) q[i] += qq[i]; // Accumulation into global q
    }
    // Note: If the outer block itself were parallelized by distributing k,
    // qq would be thread-private, and this loop would sum the thread-private
    // qq's into the global q. With serial k, qq is just a single accumulator.

    // Get number of threads (will be 1 if not compiled with OpenMP)
#if defined(_OPENMP)
    nthreads = omp_get_num_threads();
#endif     
} // End of main computational block

    // Sum elements of the global q array
    // Small loop, not the main workload
    for (i_glob = 0; i_glob <= NQ-1; i_glob++) {
        gc = gc + q[i_glob];
    }

    // Stop main timer
    timer_stop(1);
    tm = timer_read(1);

    // Verification checks - serial
    nit = 0; // nit is not used in this version of EP?
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

    // Calculate Mops
    Mops = pow(2.0, M+1)/tm/1000000.0;
	  
    // Print results
    c_print_results("EP", CLASS, M+1, 0, 0, nit, nthreads,
		  tm, Mops, 	
		  "Random numbers generated", // Description should match actual work
		  verified, NPBVERSION, COMPILETIME,
		  CS1, CS2, CS3, CS4, CS5, CS6, CS7);

    if (TIMERS_ENABLED == TRUE) {
        // Print other timers if enabled
    }

    return 0; // Added return 0 for main
}

// Helper function pow2 used in the inner loop
double pow2( double a )
{
    return a * a;
}