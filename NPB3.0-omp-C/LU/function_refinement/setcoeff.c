#include <stdio.h>
#include <math.h> // For max

// Assume these are global variables defined elsewhere
extern double dxi, deta, dzeta;
extern double tx1, tx2, tx3;
extern double ty1, ty2, ty3;
extern double tz1, tz2, tz3;
extern int ii1, ii2, ji1, ji2, ki1, ki2;
extern double dx1, dx2, dx3, dx4, dx5;
extern double dy1, dy2, dy3, dy4, dy5;
extern double dz1, dz2, dz3, dz4, dz5;
extern double dssp;
extern double ce[5][13]; // Assuming size 5x13 from assignments

// Assume nx0, ny0, nz0 are global dimensions
extern int nx0, ny0, nz0;

/*
 * Refined function to set coefficients.
 * This function calculates and assigns values to global coefficients.
 * The calculations for scalar coefficients are kept as a sequence
 * as they are fast and have minimal dependencies resolved sequentially.
 * The array initialization for 'ce' is structured into nested loops
 * to make it amenable for parallelization with OpenMP if needed,
 * by decoupling the values from the assignment process.
 */
static void setcoeff(void) {
    // Scalar coefficients calculation and assignment
    // These calculations are fast and sequential dependencies are minimal.
    // Keeping them as a sequence is efficient.
    dxi = 1.0 / ( nx0 - 1 );
    deta = 1.0 / ( ny0 - 1 );
    dzeta = 1.0 / ( nz0 - 1 );

    tx1 = 1.0 / ( dxi * dxi );
    tx2 = 1.0 / ( 2.0 * dxi );
    tx3 = 1.0 / dxi;

    ty1 = 1.0 / ( deta * deta );
    ty2 = 1.0 / ( 2.0 * deta );
    ty3 = 1.0 / deta;

    tz1 = 1.0 / ( dzeta * dzeta );
    tz2 = 1.0 / ( 2.0 * dzeta );
    tz3 = 1.0 / dzeta;

    ii1 = 1;
    ii2 = nx0 - 2;
    ji1 = 1;
    ji2 = ny0 - 3;
    ki1 = 2;
    ki2 = nz0 - 2;

    dx1 = 0.75;
    dx2 = dx1;
    dx3 = dx1;
    dx4 = dx1;
    dx5 = dx1;

    dy1 = 0.75;
    dy2 = dy1;
    dy3 = dy1;
    dy4 = dy1;
    dy5 = dy1;

    dz1 = 1.00;
    dz2 = dz1;
    dz3 = dz1;
    dz4 = dz1;
    dz5 = dz1;

    dssp = ( fmax (dx1, fmax(dy1, dz1) ) ) / 4.0; // Use fmax for double

    // Array coefficients assignment
    // These assignments are independent of each other.
    // Structuring into loops makes it easy to parallelize using OpenMP.
    // The values are decoupled from the assignment logic itself by
    // storing them in a temporary structure.
    static const double ce_values[5][13] = {
        {2.0, 0.0, 0.0, 4.0, 5.0, 3.0, 5.0e-01, 2.0e-02, 1.0e-02, 3.0e-02, 5.0e-01, 4.0e-01, 3.0e-01},
        {1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0e+00, 1.0e-02, 3.0e-02, 2.0e-02, 4.0e-01, 3.0e-01, 5.0e-01},
        {2.0, 2.0, 0.0, 0.0, 0.0, 2.0, 3.0e+00, 4.0e-02, 3.0e-02, 5.0e-02, 3.0e-01, 5.0e-01, 4.0e-01},
        {2.0, 2.0, 0.0, 0.0, 0.0, 2.0, 3.0e+00, 3.0e-02, 5.0e-02, 4.0e-02, 2.0e-01, 1.0e-01, 3.0e-01},
        {5.0, 4.0, 3.0, 2.0, 1.0e-01, 4.0e-01, 3.0e-01, 5.0e-02, 4.0e-02, 3.0e-02, 1.0e-01, 3.0e-01, 2.0e-01}
    };

    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 13; ++j) {
            ce[i][j] = ce_values[i][j];
        }
    }
}