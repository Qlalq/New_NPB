#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h> // Include OpenMP header

// Assume sprnvc, vecset, sparse are declared elsewhere and are available.
// Their internal logic is not modified here, but their thread-safety and
// reliance on scratch space require managing thread-private data.
// Dummy declarations for compilation if needed.
// Note: The original NAS SP uses 1-based indexing extensively.
// The provided snippet implies 1-based indexing for colidx, v, iv, arow, acol, aelt, rowstr.
// The refined code aims to use 1-based indexing for these arrays where they are inputs/outputs
// to match the likely expectation of sprnvc, vecset, and sparse.
// Thread-local counters and offsets will use 0-based indexing internally before converting
// to 1-based global indices for array writes.

// Dummy declarations (replace with actual headers if available)
// Assuming v, iv, colidx, a, arow, acol, aelt, rowstr are 1-based arrays
void sprnvc(int n, int nzv, double v[], int iv[], int iseed[], int nzloc[]);
void vecset(int n, double v[], int iv[], int *nzv, int i, double val);
void sparse(double a[], int colidx[], int rowstr[], int n, int arow[], int acol[], double aelt[], int firstrow, int lastrow, double v[], int iv[], int nzloc[], int nnza);


static void makea_refined(
    int n,          // Matrix size (max 1-based index is n)
    int nz,         // Max non-zeros (max 1-based index is nz)
    double a[],		// Output: dense representation (not used here, size >= n+nz+1 expected)
    int colidx[],	// Output: column indices for sparse matrix (CSR format) + scratch space. Size >= n+nz+1 expected for CSR, and >= 2*n+1 needed for scratch based on original code's use. Total size >= max(n+nz, 2*n)+1 ? Let's assume 2*n+1 is sufficient scratch, and the rest is for sparse output.
    int rowstr[],	// Output: row pointers for sparse matrix (CSR format) - use 1-based? Size >= n+1 expected.
    int nonzer,     // Number of non-zeros per generated sparse row
    int firstrow,   // First row index (1-based)
    int lastrow,    // Last row index (1-based)
    int firstcol,   // First col index (1-based)
    int lastcol,    // Last col index (1-based)
    double rcond,   // Reciprocal condition number
    double aelt[],	// Output: values (COO format initially, size >= nz+1, 1-based index)
    int arow[],		// Output: row indices (COO format initially, size >= nz+1, 1-based index)
    int acol[],		// Output: col indices (COO format initially, size >= nz+1, 1-based index)
    double v_global[], // Scratchpad/Input: vector (size >= n+1, 1-based index) - Used by sparse
    int iv_global[],   // Scratchpad/Input: vector indices (size >= n+1, 1-based index) - Used by sparse and sprnvc/vecset scratch
    double shift )  // Shift value
{
    int i, iouter;
    double ratio;
    int total_nnza_offdiag = 0;
    int total_nnza_diag = 0;
    int total_nnza;

    // 1. Initialization Loop (Zero out colidx from n+1 to 2n, using 1-based index)
    // This part of colidx is used as scratch for sprnvc/vecset.
    // This loop writes to independent locations and can be parallelized.
    #pragma omp parallel for
    for (i = 1; i <= n; i++) {
	colidx[n+i] = 0; // Use 1-based indexing as in original
    }

    ratio = pow(rcond, (1.0 / (double)n));

    int max_threads = omp_get_max_threads();

    // Thread-private counters storage (will store 0-based counts)
    // Allocate on heap to avoid potential stack overflow for large number of threads
    int* thread_offdiag_counts = (int*)calloc(max_threads, sizeof(int));
    int* thread_diag_counts = (int*)calloc(max_threads, sizeof(int));
    // Thread-private offset storage (will store 0-based starting offsets in the final COO arrays)
    int* offdiag_start_idx = (int*)calloc(max_threads, sizeof(int));
    int* diag_start_idx = (int*)calloc(max_threads, sizeof(int));

    // Check for allocation errors
    if (!thread_offdiag_counts || !thread_diag_counts || !offdiag_start_idx || !diag_start_idx) {
        fprintf(stderr, "Memory allocation failed in makea_refined\n");
        exit(1); // Exit or handle error appropriately
    }

    // --- Pass 1: Parallel Pre-counting ---
    // Each thread counts how many off-diagonal and diagonal elements it would generate.
    // Requires thread-private scratch arrays for sprnvc/vecset as they modify v and iv.
    #pragma omp parallel
    {
        int my_id = omp_get_thread_num();
        int my_nnza_offdiag = 0;
        int my_nnza_diag = 0;

        // Allocate thread-private scratch arrays (using stack for example, heap better for large n)
        // These need to be large enough for 1-based indexing up to n or 2*n.
        double v_priv[n+1];
        int iv_priv[n+1];
        int colidx_scratch_priv[2*n+1]; // Based on zeroing loop colidx[n+1..2n] and sprnvc args

        int nzv_priv; // thread-private nzv
        int ivelt, ivelt1, jcol, irow; // thread-private loop indices/variables
        double scale_priv; // thread-private scale

        // Parallelize the iouter loop for off-diagonal count
        // Static schedule ensures a predictable distribution of iterations
        #pragma omp for schedule(static) nowait // Use nowait as diagonal loop is independent
        for (iouter = 1; iouter <= n; iouter++) { // 1-based loop
            double current_size = pow(ratio, (double)(iouter - 1)); // Calculate size independently

            nzv_priv = nonzer; // nzv is constant per outer loop iteration (if nonzer is constant)

            // Call sprnvc and vecset using thread-private scratch arrays.
            // Pass addresses corresponding to 1-based indexing if sprnvc/vecset expect it.
            // Original call was &(colidx[0]), &(colidx[n]). Assuming this maps to 1-based start and 1-based start+n.
            // Pass &(colidx_scratch_priv[1]) and &(colidx_scratch_priv[n+1]) for 1-based indexing.
            // Assuming sprnvc/vecset are re-entrant and their modifications to v_priv/iv_priv/colidx_scratch_priv are thread-local.
            sprnvc(n, nzv_priv, v_priv, iv_priv, &(colidx_scratch_priv[1]), &(colidx_scratch_priv[n+1]));
            vecset(n, v_priv, iv_priv, &nzv_priv, iouter, 0.5); // nzv_priv might be modified here

            // --- Nested loops for off-diagonal count ---
            // These loops iterate over the elements produced by sprnvc/vecset for this iouter
            for (ivelt = 1; ivelt <= nzv_priv; ivelt++) { // 1-based loop
                jcol = iv_priv[ivelt];
                if (jcol >= firstcol && jcol <= lastcol) {
                    scale_priv = current_size * v_priv[ivelt];
                    for (ivelt1 = 1; ivelt1 <= nzv_priv; ivelt1++) { // 1-based loop
                        irow = iv_priv[ivelt1];
                        if (irow >= firstrow && irow <= lastrow) {
                            my_nnza_offdiag++; // Increment thread-private counter
                        }
                    }
                }
            }
        } // End parallel for iouter

        // Parallelize the i loop for diagonal count
        #pragma omp for schedule(static)
        for (i = firstrow; i <= lastrow; i++) { // 1-based loop
            if (i >= firstcol && i <= lastcol) {
                my_nnza_diag++; // Increment thread-private counter
            }
        } // End parallel for i

        // Store thread-private counts in the shared arrays
        thread_offdiag_counts[my_id] = my_nnza_offdiag;
        thread_diag_counts[my_id] = my_nnza_diag;

        // Note: In a real application, managing thread-private scratch arrays
        // dynamically using thread-local storage or master/sections might be preferred
        // over fixed-size stack allocation if N is very large.
    } // End OpenMP Parallel Region for Pass 1


    // --- Pass 2: Sequential Offset Calculation and Total Count Check ---
    // Sum thread counts to get total counts
    total_nnza_offdiag = 0;
    total_nnza_diag = 0;
    // Calculate prefix sums for 0-based offsets in the final combined COO array
    offdiag_start_idx[0] = 0; // Thread 0 starts at index 0 for off-diagonal
    diag_start_idx[0] = 0;    // Thread 0 starts at index 0 within its diagonal block

    for(int tid = 0; tid < max_threads; ++tid) {
        total_nnza_offdiag += thread_offdiag_counts[tid];
        total_nnza_diag += thread_diag_counts[tid];
        if (tid < max_threads - 1) {
            // The starting index for thread tid+1 is the starting index of thread tid + count of thread tid
            offdiag_start_idx[tid+1] = offdiag_start_idx[tid] + thread_offdiag_counts[tid];
            diag_start_idx[tid+1] = diag_start_idx[tid] + thread_diag_counts[tid];
        }
    }
    total_nnza = total_nnza_offdiag + total_nnza_diag; // Total elements = off-diagonal + diagonal

    // Check maximum size BEFORE starting the generation phase
    // This replaces the original error check inside the loop.
    if (total_nnza > nz) {
        fprintf(stderr, "Space for matrix elements exceeded in makea\n");
        fprintf(stderr, "nnza, nzmax = %d, %d\n", total_nnza, nz);
        // Original code printed iouter, which is not meaningful here.
        // Exit is necessary as the output arrays are not large enough.
        exit(1);
    }


    // --- Pass 3: Parallel Generation and Copy ---
    // Each thread generates its elements and writes them directly into the global
    // arow, acol, aelt arrays using the pre-calculated offsets.
    #pragma omp parallel
    {
        int my_id = omp_get_thread_num();
        // Base 0-based index for this thread's off-diagonal elements block
        int my_offdiag_base_0based = offdiag_start_idx[my_id];
        // Base 0-based index for this thread's diagonal elements block (starts after all offdiag)
        int my_diag_base_0based = total_nnza_offdiag + diag_start_idx[my_id];

        // Counters for elements written *within this thread's block* for current parallel loop.
        // These reset for each #pragma omp for loop but track position within the thread's assigned iterations for that loop.
        // Let's use a single counter per element type for each thread in this parallel region.
        // Initialize before the first #pragma omp for within this parallel region.
        int my_current_offdiag_count = 0;
        int my_current_diag_count = 0;


        // Allocate thread-private scratch arrays (using stack for example, heap better for large n)
        double v_priv[n+1];
        int iv_priv[n+1];
        int colidx_scratch_priv[2*n+1]; // Based on zeroing loop colidx[n+1..2n] and sprnvc args

        int nzv_priv;
        int ivelt, ivelt1, jcol, irow;
        double scale_priv;


        // Parallel loop for off-diagonal generation and writing
        // Use same schedule as in Pass 1 to ensure same iteration distribution
        #pragma omp for schedule(static) nowait
        for (iouter = 1; iouter <= n; iouter++) { // 1-based loop
            double current_size = pow(ratio, (double)(iouter - 1));

            nzv_priv = nonzer;
            // Call sprnvc and vecset using thread-private scratch
            sprnvc(n, nzv_priv, v_priv, iv_priv, &(colidx_scratch_priv[1]), &(colidx_scratch_priv[n+1]));
            vecset(n, v_priv, iv_priv, &nzv_priv, iouter, 0.5);

            for (ivelt = 1; ivelt <= nzv_priv; ivelt++) { // 1-based loop
                jcol = iv_priv[ivelt];
                if (jcol >= firstcol && jcol <= lastcol) {
                    scale_priv = current_size * v_priv[ivelt];
                    for (ivelt1 = 1; ivelt1 <= nzv_priv; ivelt1++) { // 1-based loop
                        irow = iv_priv[ivelt1];
                        if (irow >= firstrow && irow <= lastrow) {
                            // Calculate the 0-based write index in the combined COO arrays
                            int write_idx_0based = my_offdiag_base_0based + my_current_offdiag_count;

                            // Write using 1-based indexing in the global arow, acol, aelt arrays
                            // (which are assumed to be 1-based, size nz+1)
                            arow[write_idx_0based + 1] = irow; // Store 1-based row index
                            acol[write_idx_0based + 1] = jcol; // Store 1-based col index
                            aelt[write_idx_0based + 1] = v_priv[ivelt1] * scale_priv;

                            my_current_offdiag_count++; // Increment thread-local counter for this block
                        }
                    }
                }
            }
        } // End parallel for iouter

        // Parallel loop for diagonal generation and writing
        // Use same schedule as in Pass 1
        #pragma omp for schedule(static) nowait
        for (i = firstrow; i <= lastrow; i++) { // 1-based loop
            if (i >= firstcol && i <= lastcol) {
                // Calculate the 0-based write index in the combined COO arrays
                int write_idx_0based = my_diag_base_0based + my_current_diag_count;

                // Write using 1-based indexing in the global arow, acol, aelt arrays
                arow[write_idx_0based + 1] = i; // Store 1-based row index
                acol[write_idx_0based + 1] = i; // Store 1-based col index
                aelt[write_idx_0based + 1] = rcond - shift;

                my_current_diag_count++; // Increment thread-local counter for this block
            }
        } // End parallel for i

        // Note: Thread-private scratch arrays allocated on stack are freed automatically here.
        // If using dynamic allocation, free them before the end of the parallel region.

    } // End OpenMP Parallel Region for Pass 3


    // --- Pass 4: Sequential Finalization ---
    // Call sparse to convert COO (arow, acol, aelt) to CSR (a, colidx, rowstr).
    // sparse expects arow, acol, aelt to be 1-based arrays, size total_nnza.
    // (although declared size nz+1, only total_nnza elements are valid, indices 1..total_nnza).
    // sparse also uses v_global and iv_global as workspace (likely size n+1, 1-based).
    // It takes nzloc as &(iv_global[n+1]) based on typical SP usage.
    // It takes nnza as the number of elements (total_nnza).
    sparse(a, colidx, rowstr, n, arow, acol, aelt,
	   firstrow, lastrow, v_global, &(iv_global[1]), &(iv_global[n+1]), total_nnza);


    // Free allocated memory for counts and offsets
    free(thread_offdiag_counts);
    free(thread_diag_counts);
    free(offdiag_start_idx);
    free(diag_start_idx);
}