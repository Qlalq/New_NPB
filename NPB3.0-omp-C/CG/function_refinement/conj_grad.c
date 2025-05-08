#include <math.h> // Needed for sqrt

// Assume naa, lastcol, firstcol, lastrow, firstrow are defined elsewhere
// (e.g., as global variables or preprocessor macros)
extern int naa;
extern int lastcol;
extern int firstcol;
extern int lastrow;
extern int firstrow;


static void conj_grad_refined (
    int colidx[],	
    int rowstr[],	
    double x[],		
    double z[],		
    double a[],		
    double p[],		
    double q[],		
    double r[],		
    //double w[], // w is commented out in original, omit
    double *rnorm )
{
    static int callcount = 0; // static shared variable, increment needs care in parallel
    double d;
    double sum; // Used for local sums in A*p, A*z and global sum for rnorm
    double rho, rho0, alpha, beta;
    int i, j, k;
    int cgit, cgitmax = 25;

    // ---------------------------------------------------------------------
    // 1. Initialization: z = 0, r = x, p = x, q = 0
    //    This section consists of element-wise vector assignments.
    //    Each iteration is independent.
    // ---------------------------------------------------------------------
    for (j = 1; j <= naa+1; j++) {
        q[j] = 0.0;
    }

    for (j = 1; j <= lastcol-firstcol+1; j++) {
        z[j] = 0.0;
        r[j] = x[j];
        p[j] = r[j];
        //w[j] = 0.0; // w is commented out
    }

    // ---------------------------------------------------------------------
    // 2. Calculate initial rho = r.r (dot product)
    //    This is a reduction operation.
    // ---------------------------------------------------------------------
    rho = 0.0; // Reduction variable initialization
    for (j = 1; j <= lastcol-firstcol+1; j++) {
        rho += r[j] * r[j]; // Summation loop
    }


    // ---------------------------------------------------------------------
    // Main Conjugate Gradient iteration loop
    // ---------------------------------------------------------------------
    for (cgit = 1; cgit <= cgitmax; cgit++) {

        rho0 = rho; // Save rho from previous iteration (scalar assignment)
        d = 0.0; // Initialize d for the p.q reduction

        // ---------------------------------------------------------------------
        // 3. Calculate q = A*p (sparse matrix-vector multiplication)
        //    The outer loop over rows (j) can be parallelized.
        //    The inner loop computes a local sum for each row.
        // ---------------------------------------------------------------------
        for (j = 1; j <= lastrow-firstrow+1; j++) {
            sum = 0.0; // Local sum for each row j - should be thread private
            for (k = rowstr[j]; k < rowstr[j+1]; k++) {
                sum += a[k] * p[colidx[k]]; // Inner loop dependency on local sum
            }
            q[j] = sum; // Independent write for each j
        }

        // ---------------------------------------------------------------------
        // 4. Calculate d = p.q (dot product)
        //    This is a reduction operation.
        // ---------------------------------------------------------------------
        d = 0.0; // Reduction variable initialization
        for (j = 1; j <= lastcol-firstcol+1; j++) {
            d += p[j] * q[j]; // Summation loop
        }

        // ---------------------------------------------------------------------
        // 5. Calculate alpha
        //    Scalar assignment, depends on results of reductions.
        // ---------------------------------------------------------------------
        alpha = rho0 / d;

        // ---------------------------------------------------------------------
        // 6. Update z and r vectors: z = z + alpha*p, r = r - alpha*q
        //    These are element-wise vector operations. Each iteration is independent.
        // ---------------------------------------------------------------------
        for (j = 1; j <= lastcol-firstcol+1; j++) {
            z[j] = z[j] + alpha * p[j];
            r[j] = r[j] - alpha * q[j];
        }

        // ---------------------------------------------------------------------
        // 7. Calculate new rho = r.r (dot product)
        //    This is a reduction operation on the *updated* r vector.
        // ---------------------------------------------------------------------
        rho = 0.0; // Reduction variable initialization for the next iteration's rho0
        for (j = 1; j <= lastcol-firstcol+1; j++) {
            rho += r[j] * r[j]; // Summation loop
        }

        // ---------------------------------------------------------------------
        // 8. Calculate beta
        //    Scalar assignment, depends on results of reductions.
        // ---------------------------------------------------------------------
        beta = rho / rho0;

        // ---------------------------------------------------------------------
        // 9. Update p vector: p = r + beta*p
        //    This is an element-wise vector operation. Each iteration is independent.
        // ---------------------------------------------------------------------
        for (j = 1; j <= lastcol-firstcol+1; j++) {
            p[j] = r[j] + beta * p[j];
        }

        // ---------------------------------------------------------------------
        // Call count update (scalar, requires atomic/critical in parallel)
        // ---------------------------------------------------------------------
        callcount++;

    } // End of main CG iteration loop

    // ---------------------------------------------------------------------
    // Post-iteration calculations
    // ---------------------------------------------------------------------

    // ---------------------------------------------------------------------
    // 10. Calculate final r = A*z (sparse matrix-vector multiplication)
    //     Similar structure to step 3.
    // ---------------------------------------------------------------------
    for (j = 1; j <= lastrow-firstrow+1; j++) {
        sum = 0.0; // Local sum for each row j - should be thread private
        // Note: Original code used d here as local sum, switching to sum for consistency with step 3
        for (k = rowstr[j]; k <= rowstr[j+1]-1; k++) { // Original loop bounds
            sum += a[k] * z[colidx[k]];
        }
        r[j] = sum; // Independent write for each j
    }

    // ---------------------------------------------------------------------
    // 11. Calculate sum = (x-r).(x-r) (dot product of difference)
    //     This is a reduction operation.
    // ---------------------------------------------------------------------
    sum = 0.0; // Reduction variable initialization
    for (j = 1; j <= lastcol-firstcol+1; j++) {
        d = x[j] - r[j]; // Local temporary for the difference
        sum += d * d; // Summation loop (d is used locally here, no conflict with d from step 4)
    }

    // ---------------------------------------------------------------------
    // 12. Calculate final rnorm
    //     Scalar assignment, depends on the final reduction result.
    // ---------------------------------------------------------------------
    (*rnorm) = sqrt(sum);

}