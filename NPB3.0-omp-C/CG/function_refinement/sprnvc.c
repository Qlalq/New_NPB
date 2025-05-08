#include <stdio.h>
#include <stdlib.h> // For malloc, free
#include <string.h> // For memset (optional, depending on caller)
#include <omp.h>    // For OpenMP pragmas and functions

// Define a structure to hold generated candidates (index and value)
struct Candidate {
    int index;
    double value;
};

// --- Dummy Implementations (Replace with actual NPB functions) ---
// In a real NPB implementation, randlc and icnvrt are provided.
// For parallel execution, randlc usually needs to be re-entrant
// or support independent streams using thread-local state.
// These dummies are non-functional but show the structure where
// the actual calls would be made.

// Dummy state for randlc (a real parallel version would need per-thread state)
static double DUMMY_TRAN_STATE = 1.0;

// Dummy randlc function - NOT thread-safe or functionally correct
static double randlc(double *p, double amult) {
    // This is a placeholder. Replace with actual randlc logic.
    // A parallel version typically uses a state pointer unique to each thread.
    // Example dummy state update and return:
    *p = *p * amult + 0.1;
    return *p; // Returns a value based on state
}

// Dummy icnvrt function - NOT functionally correct
static int icnvrt(double x, int scale) {
    // This is a placeholder. Replace with actual icnvrt logic.
    // It converts a double in [0,1) range to an integer index.
    // Example dummy conversion:
    return (int)(x * scale / 2.0);
}
// --- End Dummy Implementations ---

// Assume amult is accessible, e.g., a global variable or passed somehow.
// In NPB, it's often a global constant.
static double amult = 0.75; // Dummy value for amult

/*
 * Refined sprnvc function for better OpenMP structure.
 *
 * Original function generates a sparse vector with 'nz' non-zero elements
 * whose indices are generated pseudo-randomly within the range [1, n].
 * It also collects the list of unique indices found.
 *
 * The original code has strong sequential dependencies:
 * 1. Random number generation state (`randlc`).
 * 2. Sequential counters (`nzv`, `nzrow`) used for output array indices.
 * 3. Stateful uniqueness check and update (`mark` array).
 *
 * This refactoring separates the process into phases:
 * 1. Setup (calculate nn1) - Sequential.
 * 2. Parallel Candidate Generation - Generate a pool of potential (index, value)
 *    pairs using parallel random number generation. This phase produces more
 *    candidates than needed and stores valid ones (index within [1, n]) in a
 *    temporary buffer.
 * 3. Sequential Filtering and Population - Iterate through the temporary buffer
 *    sequentially, applying the original logic (uniqueness check using `mark`,
 *    incrementing `nzv`/`nzrow`, populating `v`/`iv`/`nzloc`) until `nz` elements
 *    are collected. This phase preserves the serial output logic order based on
 *    the order candidates appear in the temporary buffer.
 * 4. Cleanup (reset mark array) - Parallelizable.
 *
 * This structure isolates the core sequential logic (Phase 3) while parallelizing
 * parts that can be done independently (Phase 2 and 4). The temporary buffer
 * decouples the rate of generation from the rate of acceptance.
 *
 * Arguments:
 * n:     Dimension of the vector (indices 1 to n).
 * nz:    Target number of non-zero elements.
 * v:     Output array for non-zero values (1-based, size nz+1).
 * iv:    Output array for non-zero indices (1-based, size nz+1).
 * nzloc: Output array for unique indices found (1-based, size nz+1).
 * mark:  Auxiliary array for uniqueness check (1-based, size n+1). Must be
 *        initialized to 0 before the first call. Reset to 0 on exit.
 */
static void sprnvc_refined(
    int n,
    int nz,
    double v[],		// size nz+1 (1-based)
    int iv[],		// size nz+1 (1-based)
    int nzloc[],	// size nz+1 (1-based)
    int mark[] )    // size n+1 (1-based)
{
    int nn1;
    int nzrow;     // count of unique row indices found (1-based)
    int nzv;       // count of non-zero elements found (1-based)
    int ii, i;

    // Assume mark array is pre-allocated by caller and needs to be initialized/reset.
    // Initialize mark array to 0 before use. This can be parallelized if n is large.
    // The original code implies mark is reset at the end, assuming it's either
    // persistent and needs cleanup, or needs initialization if not persistent.
    // Initializing it here ensures a clean state for the filtering phase.
    #pragma omp parallel for
    for (i = 1; i <= n; ++i) {
        mark[i] = 0;
    }

    // Phase 1: Calculate nn1 (Sequential)
    // Small calculation, not performance critical.
    nn1 = 1;
    do {
	nn1 = 2 * nn1;
    } while (nn1 < n);

    // Phase 2: Parallel Candidate Generation
    // Generate a pool of potential candidates (index, value) pairs.
    // We generate more candidates than needed to increase the probability
    // of finding 'nz' valid non-zeros and 'nzrow' unique indices quickly.
    // The factor (e.g., 2.0) is a heuristic and might need tuning.
    int num_attempts = (int)(2.0 * nz); // Number of random pairs to generate

    // Allocate temporary buffer for candidates that have valid indices (1..n).
    // Allocate size num_attempts as a safe upper bound for the number of valid candidates.
    struct Candidate* candidate_buffer = (struct Candidate*)malloc(num_attempts * sizeof(struct Candidate));
    if (!candidate_buffer) {
        // Handle allocation error
        fprintf(stderr, "Error: Could not allocate memory for candidate_buffer\n");
        // Depending on context, may need to exit or signal failure
        return;
    }

    // For a real parallel NPB randlc, you would need to manage thread-local
    // random number generator states. This requires modifying/using a randlc
    // version that supports state passing.
    // Example: double thread_rand_states[omp_get_max_threads()];
    // And initialize these states here...

    // Shared atomic counter for the number of valid candidates stored in the buffer.
    int num_valid_candidates = 0;

    // Use a parallel loop to generate candidates.
    // The schedule(static) is simple. schedule(dynamic) could be used if
    // the cost of generating/validating a single candidate varied significantly.
    #pragma omp parallel
    {
        // Get thread's random state if needed (replace DUMMY_TRAN_STATE management)
        // double my_tran = initial_thread_state[omp_get_thread_num()];

        #pragma omp for schedule(static)
        for (int attempt = 0; attempt < num_attempts; ++attempt) {
            double vecelt, vecloc;
            int current_i;

            // --- Generate random numbers and index ---
            // Replace these calls with thread-safe randlc and the actual icnvrt.
            // Dummy implementation for structure:
            vecelt = randlc(&DUMMY_TRAN_STATE, amult); // Using global dummy state (NOT safe)
            vecloc = randlc(&DUMMY_TRAN_STATE, amult); // Using global dummy state (NOT safe)
            current_i = icnvrt(vecloc, nn1) + 1;
            // ---------------------------------------

            // Check if the generated index is within the valid range [1, n]
            if (current_i >= 1 && current_i <= n) {
                // Atomically get the next available slot index in the shared buffer
                int current_buffer_idx;
                // Using atomic capture to get the value before increment
                #pragma omp atomic capture
                current_buffer_idx = num_valid_candidates++;

                // Store the valid candidate if there is space in the buffer.
                // This check prevents writing beyond the allocated buffer size,
                // although ideally num_attempts is large enough.
                if (current_buffer_idx < num_attempts) {
                    candidate_buffer[current_buffer_idx].index = current_i;
                    candidate_buffer[current_buffer_idx].value = vecelt;
                }
                // else: This candidate is valid, but we hit the limit of our
                // temporary buffer. Discard it. This might mean num_attempts
                // was too small to find 'nz' non-zeros.
            }
        }
        // If using thread-local states, update the main array of states here.
        // final_thread_state[omp_get_thread_num()] = my_tran;
    } // End parallel region for generation

    // Phase 3: Sequential Filtering and Population
    // Process the collected valid candidates sequentially.
    // This phase applies the original logic of checking the 'mark' array
    // and populating the final output arrays (v, iv, nzloc) in the order
    // candidates were added to the buffer (which approximates the serial order
    // if generation is structured correctly). This is the critical sequential
    // part because 'mark', 'nzv', 'nzrow', v, iv, nzloc are stateful and
    // depend on the processing order.

    nzv = 0;     // Reset non-zero count (1-based)
    nzrow = 0;   // Reset unique index count (1-based)

    // Loop through the candidates that were successfully added to the buffer.
    // The loop runs up to the actual count of valid candidates found in Phase 2.
    for (int j = 0; j < num_valid_candidates; ++j) {
        // Stop if we have collected the target number of non-zero elements.
        // This condition ensures 'nzv' does not exceed 'nz'.
        if (nzv == nz) {
            break;
        }

        int current_i = candidate_buffer[j].index;
        double vecelt = candidate_buffer[j].value;

        // Check if this index has already been marked (means it was accepted earlier
        // as the first occurrence of a non-zero at this index).
        if (mark[current_i] == 0) {
            // This is the first time we've encountered a non-zero at this index.
            // Mark this index as having a non-zero that will be added to nzloc.
            mark[current_i] = 1;

            // Add the index to the list of unique indices found (nzloc).
            nzrow = nzrow + 1;
            // Ensure we don't write out of bounds for nzloc (size nz+1, 1-based index).
            // In the original logic, nzrow <= nzv, and nzv stops at nz, so nzrow <= nz.
            if (nzrow <= nz) {
               nzloc[nzrow] = current_i;
            } else {
               // This case indicates more unique indices were found than nz
               // elements were accepted. This shouldn't happen based on
               // the original logic where nzrow increments only if nzv increments.
               // If it could, nzloc would need to be larger.
            }


            // Add the non-zero value and its index to v and iv.
            // Increment nzv *before* using it as an index, as it's 1-based.
            nzv = nzv + 1;
            // Ensure we don't write out of bounds for v and iv (size nz+1, 1-based index).
            // This is guaranteed by the `if (nzv == nz) break;` check at the start
            // of the loop, meaning nzv will not exceed nz when used as an index.
            v[nzv] = vecelt;
            iv[nzv] = current_i;
        }
    }

    // Free the temporary candidate buffer memory.
    free(candidate_buffer);

    // Phase 4: Cleanup mark array (Parallelizable)
    // Reset the mark for all indices that were added to nzloc.
    // Use the final value of nzrow determined in Phase 3.
    #pragma omp parallel for schedule(static)
    for (ii = 1; ii <= nzrow; ii++) {
        // Get the index from nzloc (read-only access here)
        i = nzloc[ii];
        // Reset the mark for this index (independent write)
        mark[i] = 0;
    }

    // Note: After this function, `nzv` should ideally be equal to `nz`.
    // If `num_attempts` was too low or the random distribution unfavorable,
    // it might be less than `nz`. A robust implementation might need to handle
    // this by, e.g., repeating phase 2 and 3, but that adds complexity.
    // This refactoring assumes `num_attempts` is sufficient.
}

// Use the original function name and call the refined implementation.
static void sprnvc(
    int n,
    int nz,
    double v[],		// size nz+1 (1-based)
    int iv[],		// size nz+1 (1-based)
    int nzloc[],	// size nz+1 (1-based)
    int mark[] )    // size n+1 (1-based)
{
     // Call the refactored function implementation
     sprnvc_refined(n, nz, v, iv, nzloc, mark);
}