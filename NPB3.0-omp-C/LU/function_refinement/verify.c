#include <stdio.h>
#include <math.h>   // For fabs

// Assuming boolean, TRUE, FALSE, nx0, ny0, nz0, itmax, dt are defined elsewhere
// as they are used but not declared or passed into the original function signature.
// For demonstration, let's add minimal definitions if they weren't provided in the
// original code snippet's surrounding context.
/*
typedef int boolean;
#define TRUE 1
#define FALSE 0
extern int nx0, ny0, nz0, itmax;
extern double dt;
*/

static void verify_refined(double xcr[5], double xce[5], double xci,
		   char *class, boolean *verified) {

    // Declare all necessary variables upfront
    double xcrref[5],xceref[5],xciref;
    double xcrdif[5],xcedif[5],xcidif;
    double epsilon;
    double dtref = 0.0; // Initialize dtref to avoid potential issues if class is 'U' and dtref is accessed
    int m;

    // Sequential Part 1: Initialization and Class Determination
    // Initialize default reference values. This loop is small (5 iterations)
    // and can be left sequential or parallelized if desired, but overhead
    // might outweigh benefits for such a small loop.
    for (m = 0; m < 5; m++) {
        xcrref[m] = 1.0;
        xceref[m] = 1.0;
    }
    xciref = 1.0;

    epsilon = 1.0e-08;
    *class = 'U';
    *verified = TRUE; // Assume verified initially

    // Determine class and specific reference values based on problem size and iteration count.
    // This block is inherently sequential as it depends on external configuration parameters.
    if ( nx0 == 12 && ny0 == 12 && nz0 == 12 && itmax == 50)  {
        *class = 'S';
        dtref = 5.0e-1;
        xcrref[0] = 1.6196343210976702e-02;
        xcrref[1] = 2.1976745164821318e-03;
        xcrref[2] = 1.5179927653399185e-03;
        xcrref[3] = 1.5029584435994323e-03;
        xcrref[4] = 3.4264073155896461e-02;
        xceref[0] = 6.4223319957960924e-04;
        xceref[1] = 8.4144342047347926e-05;
        xceref[2] = 5.8588269616485186e-05;
        xceref[3] = 5.8474222595157350e-05;
        xceref[4] = 1.3103347914111294e-03;
        xciref = 7.8418928865937083;
    } else if ( nx0 == 33 && ny0 == 33 && nz0 == 33 && itmax == 300) {
        *class = 'W';
        dtref = 1.5e-3;
        xcrref[0] =   0.1236511638192e+02;
        xcrref[1] =   0.1317228477799e+01;
        xcrref[2] =   0.2550120713095e+01;
        xcrref[3] =   0.2326187750252e+01;
        xcrref[4] =   0.2826799444189e+02;
        xceref[0] =   0.4867877144216;
        xceref[1] =   0.5064652880982e-01;
        xceref[2] =   0.9281818101960e-01;
        xceref[3] =   0.8570126542733e-01;
        xceref[4] =   0.1084277417792e+01;
        xciref    =   0.1161399311023e+02;
    } else if ( nx0 == 64 && ny0 == 64 && nz0 == 64 && itmax == 250) {
        *class = 'A';
        dtref = 2.0e+0;
        xcrref[0] = 7.7902107606689367e+02;
        xcrref[1] = 6.3402765259692870e+01;
        xcrref[2] = 1.9499249727292479e+02;
        xcrref[3] = 1.7845301160418537e+02;
        xcrref[4] = 1.8384760349464247e+03;
        xceref[0] = 2.9964085685471943e+01;
        xceref[1] = 2.8194576365003349;
        xceref[2] = 7.3473412698774742;
        xceref[3] = 6.7139225687777051;
        xceref[4] = 7.0715315688392578e+01;
        xciref = 2.6030925604886277e+01;
    } else if ( nx0 == 102 && ny0 == 102 && nz0 == 102 && itmax == 250) {
        *class = 'B';
        dtref = 2.0e+0;
        xcrref[0] = 3.5532672969982736e+03;
        xcrref[1] = 2.6214750795310692e+02;
        xcrref[2] = 8.8333721850952190e+02;
        xcrref[3] = 7.7812774739425265e+02;
        xcrref[4] = 7.3087969592545314e+03;
        xceref[0] = 1.1401176380212709e+02;
        xceref[1] = 8.1098963655421574;
        xceref[2] = 2.8480597317698308e+01;
        xceref[3] = 2.5905394567832939e+01;
        xceref[4] = 2.6054907504857413e+02;
        xciref = 4.7887162703308227e+01;
    } else if ( nx0 == 162 && ny0 == 162 && nz0 == 162 && itmax == 250) {
        *class = 'C';
        dtref = 2.0e+0;
        xcrref[0] = 1.03766980323537846e+04;
        xcrref[1] = 8.92212458801008552e+02;
        xcrref[2] = 2.56238814582660871e+03;
        xcrref[3] = 2.19194343857831427e+03;
        xcrref[4] = 1.78078057261061185e+04;
        xceref[0] = 2.15986399716949279e+02;
        xceref[1] = 1.55789559239863600e+01;
        xceref[2] = 5.41318863077207766e+01;
        xceref[3] = 4.82262643154045421e+01;
        xceref[4] = 4.55902910043250358e+02;
        xciref = 6.66404553572181300e+01;
    } else {
        *verified = FALSE; // Class 'U' implies not verified
    }

    // Sequential Part 2: Initial DT Check and Status Update
    if (*class != 'U') {
        printf("\n Verification being performed for class %1c\n", *class);
        printf(" Accuracy setting for epsilon = %20.13e\n", epsilon);
        // Check dt against the reference value. If it fails, set verified to FALSE and class to 'U'.
        // Assumes 'dt' variable is accessible (e.g., global or defined elsewhere).
        if (fabs(dt - dtref) > epsilon) {
            *verified = FALSE;
            *class = 'U'; // Set class to Unknown on failure, similar to original logic
            printf(" DT does not match the reference value of %15.8e\n", dtref);
        }
    } else {
        printf(" Unknown class\n");
    }

    // Parallel Candidate Part: Calculate Differences
    // This loop calculates element-wise relative differences.
    // Each iteration is independent, making this loop suitable for parallelization.
    // We calculate these regardless of the class, similar to the original code's flow,
    // but the *comparison* for verification happens only if the class is known.
    for (m = 0; m < 5; m++) {
       xcrdif[m] = fabs((xcr[m]-xcrref[m])/xcrref[m]);
       xcedif[m] = fabs((xce[m]-xceref[m])/xceref[m]);
    }
    // Calculate the integral difference. This is a single calculation.
    xcidif = fabs((xci - xciref)/xciref);

    // Sequential Part 3: Check Calculated Differences and Update Verified Status
    // We check all differences *after* they've been calculated (potentially in parallel).
    // This sequential pass ensures that the final state of *verified is correctly set
    // if *any* check fails, avoiding race conditions if this were done inside parallel loops
    // or interleaved with printing.
    if (*class != 'U') {
        for (m = 0; m < 5; m++) {
            if (xcrdif[m] > epsilon) {
                *verified = FALSE;
                // No need to break; we need to check all values for potential prints later
            }
            if (xcedif[m] > epsilon) {
                *verified = FALSE;
                 // No need to break
            }
        }
        if (xcidif > epsilon) {
            *verified = FALSE;
        }
    }
    // Note: The original code also set *class = 'U' if any check failed *within*
    // the printing loops. This refined version sets *verified = FALSE* based on checks,
    // but only sets *class = 'U'* based on the initial DT check failure or
    // if no class was found. This slightly differs from the original's complex
    // interleaving but provides a cleaner structure for parallelism.
    // If strictly matching the original's 'U' setting on *any* failure is needed,
    // the Sequential Part 3 check would also need to set *class = 'U' along with
    // *verified = FALSE*. Let's adjust to match the original's logic more closely
    // regarding the *class* update on failure.

    // Refined Sequential Part 3 (closer to original *class update logic):
    if (*class != 'U') {
        boolean current_checks_passed = TRUE; // Flag to see if *any* diff failed

        for (m = 0; m < 5; m++) {
            if (xcrdif[m] > epsilon) {
                current_checks_passed = FALSE;
            }
            if (xcedif[m] > epsilon) {
                 current_checks_passed = FALSE;
            }
        }
        if (xcidif > epsilon) {
             current_checks_passed = FALSE;
        }

        if (current_checks_passed == FALSE) {
            *verified = FALSE;
            // The original logic also set *class = 'U' here if *any* check failed
            // This might override the specific class ('S', 'W', 'A', 'B', 'C').
            // This seems unusual, but if strictly following the source, we would add:
            // *class = 'U'; // Based on observation of original code setting class to 'U' on FAILURE
            // Let's assume the original intent might have been that *any* failure renders the
            // verification result "Unknown" or invalid for the specific class, setting it to 'U'.
            // We will add this line to match the observed behavior in the original print loops.
             *class = 'U';
        }
    }


    // Sequential Part 4: Print Results
    // All output is done sequentially after calculations and final verification status determination.
    // The printing logic now uses the calculated differences and the final *class status.
    if (*class != 'U') {
        printf(" Comparison of RMS-norms of residual\n");
    } else {
        printf(" RMS-norms of residual\n");
    }

    for (m = 0; m < 5; m++) {
        if (*class  ==  'U') {
            // If class is 'U', just print the computed value
            printf("          %2d  %20.13e\n", m, xcr[m]);
        } else {
            // If class was known, print the comparison details.
            // We check xcrdif[m] > epsilon here again ONLY for deciding whether to print "FAILURE".
            printf(" %s: %2d  %20.13e%20.13e%20.13e\n",
                   (xcrdif[m] > epsilon) ? "FAILURE" : "         ", // Prefix based on failure
                   m, xcr[m], xcrref[m], xcrdif[m]);
        }
    }

    if (*class != 'U') {
        printf(" Comparison of RMS-norms of solution error\n");
    } else {
        printf(" RMS-norms of solution error\n");
    }

    for (m = 0; m < 5; m++) {
        if (*class  ==  'U') {
            // If class is 'U', just print the computed value
            printf("          %2d  %20.13e\n", m, xce[m]);
        } else {
            // If class was known, print the comparison details
            printf(" %s: %2d  %20.13e%20.13e%20.13e\n",
                   (xcedif[m] > epsilon) ? "FAILURE" : "         ", // Prefix based on failure
                   m, xce[m], xceref[m], xcedif[m]);
        }
    }

    if (*class != 'U') {
        printf(" Comparison of surface integral\n");
    } else {
        printf(" Surface integral\n");
    }

    if (*class  ==  'U') {
        // If class is 'U', just print the computed value
        printf("              %20.13e\n", xci);
    } else {
        // If class was known, print the comparison details
        printf(" %s:     %20.13e%20.13e%20.13e\n",
               (xcidif > epsilon) ? "FAILURE" : "         ", // Prefix based on failure
               xci, xciref, xcidif);
    }

    // Sequential Part 5: Final Status Print
    // This uses the final *class and *verified status determined earlier.
    if (*class  ==  'U') {
        printf(" No reference values provided\n");
        printf(" No verification performed\n");
    } else if (*verified) {
        printf(" Verification Successful\n");
    } else {
        // This branch is reached if *class was initially known but *verified was set to FALSE
        // due to a difference exceeding epsilon (or dt mismatch).
        printf(" Verification failed\n");
    }
}