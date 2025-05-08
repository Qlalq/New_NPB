#include <stdio.h>
#include <math.h>   // For fabs
// Assume boolean, TRUE, FALSE are defined elsewhere
// e.g., #define boolean int
//       #define TRUE 1
//       #define FALSE 0
// Assume get_real, get_imag are defined elsewhere and take sums[i]
// Assume sums is globally or externally accessible and indexed by i

static void verify (int d1, int d2, int d3, int nt,
		    boolean *verified, char *class) {
    int i;
    double err;
    double epsilon = 1.0e-12;
    
    // Verification data arrays are constant and read-only, posing no parallelism issues
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
    
    // A local flag to track if any error occurred within a specific verification block.
    // This flag can be the target of an OpenMP reduction (e.g., boolean OR).
    boolean block_verification_failed;

    // Initialize the overall verification result flag.
    *verified = TRUE;
    // Initialize the class flag.
    *class = 'U';

    // The if-else if structure is sequential and selects the verification case.
    if (d1 == 64 && d2 == 64 && d3 == 64 && nt == 6) {
        *class = 'S';
        // Initialize the local flag for this block.
        block_verification_failed = FALSE;
        // The loop iterations are independent. Each iteration calculates and checks
        // errors for sums[i] against fixed data.
        // The dependency on *verified and the 'break' are removed.
        // An OpenMP parallel for directive with a reduction on block_verification_failed
        // can be applied here to parallelize the loop.
        for (i = 1; i <= nt; i++) {
            // Check real part
            err = (get_real(sums[i]) - vdata_real_s[i]) / vdata_real_s[i];
            if (fabs(err) > epsilon) {
                block_verification_failed = TRUE; // Set local flag, do not break
            }
            // Check imaginary part
            err = (get_imag(sums[i]) - vdata_imag_s[i]) / vdata_imag_s[i];
            if (fabs(err) > epsilon) {
                block_verification_failed = TRUE; // Set local flag, do not break
            }
        }
        // After the loop (and potential reduction), update the overall result.
        if (block_verification_failed) {
             *verified = FALSE;
        }
    } else if (d1 == 128 && d2 == 128 && d3 == 32 && nt == 6) {
        *class = 'W';
        block_verification_failed = FALSE;
        // This loop is independent per iteration.
        // OpenMP #pragma omp parallel for reduction(||:block_verification_failed) private(i, err)
        for (i = 1; i <= nt; i++) {
            err = (get_real(sums[i]) - vdata_real_w[i]) / vdata_real_w[i];
            if (fabs(err) > epsilon) {
                block_verification_failed = TRUE;
            }
            err = (get_imag(sums[i]) - vdata_imag_w[i]) / vdata_imag_w[i];
            if (fabs(err) > epsilon) {
                block_verification_failed = TRUE;
            }
        }
         if (block_verification_failed) {
             *verified = FALSE;
        }
    } else if (d1 == 256 && d2 == 256 && d3 == 128 && nt == 6) {
        *class = 'A';
        block_verification_failed = FALSE;
        // This loop is independent per iteration.
        // OpenMP #pragma omp parallel for reduction(||:block_verification_failed) private(i, err)
        for (i = 1; i <= nt; i++) {
            err = (get_real(sums[i]) - vdata_real_a[i]) / vdata_real_a[i];
            if (fabs(err) > epsilon) {
                block_verification_failed = TRUE;
            }
            err = (get_imag(sums[i]) - vdata_imag_a[i]) / vdata_imag_a[i];
            if (fabs(err) > epsilon) {
                block_verification_failed = TRUE;
            }
        }
         if (block_verification_failed) {
             *verified = FALSE;
        }
    } else if (d1 == 512 && d2 == 256 && d3 == 256 && nt == 20) {
        *class = 'B';
        block_verification_failed = FALSE;
        // This loop is independent per iteration.
        // OpenMP #pragma omp parallel for reduction(||:block_verification_failed) private(i, err)
        for (i = 1; i <= nt; i++) {
            err = (get_real(sums[i]) - vdata_real_b[i]) / vdata_real_b[i];
            if (fabs(err) > epsilon) {
                block_verification_failed = TRUE;
            }
            err = (get_imag(sums[i]) - vdata_imag_b[i]) / vdata_imag_b[i];
            if (fabs(err) > epsilon) {
                block_verification_failed = TRUE;
            }
        }
         if (block_verification_failed) {
             *verified = FALSE;
        }
    } else if (d1 == 512 && d2 == 512 && d3 == 512 && nt == 20) {
        *class = 'C';
        block_verification_failed = FALSE;
        // This loop is independent per iteration.
        // OpenMP #pragma omp parallel for reduction(||:block_verification_failed) private(i, err)
        for (i = 1; i <= nt; i++) {
            err = (get_real(sums[i]) - vdata_real_c[i]) / vdata_real_c[i];
            if (fabs(err) > epsilon) {
                block_verification_failed = TRUE;
            }
            err = (get_imag(sums[i]) - vdata_imag_c[i]) / vdata_imag_c[i];
            if (fabs(err) > epsilon) {
                block_verification_failed = TRUE;
            }
        }
         if (block_verification_failed) {
             *verified = FALSE;
        }
    }

    // These print statements should remain after the verification logic is complete.
    if (*class != 'U') {
	printf("Result verification successful\n");
    } else {
	printf("Result verification failed\n");
    }
    printf("class = %1c\n", *class);
}