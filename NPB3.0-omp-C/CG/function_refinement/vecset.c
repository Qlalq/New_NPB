#include <stdio.h> // Assuming standard headers might be needed or were in the original context
// Assuming definition for boolean and FALSE exists, adding them here for completeness
typedef int boolean;
#define FALSE 0
#define TRUE 1

static void vecset_refined(
    int n,
    double v[],
    int iv[],
    int *nzv,
    int i,
    double val)
{
    // Objective: Refine vecset for better OpenMP parallelization potential.
    // Original logic: Search for index 'i' in iv. If found, update v[k]. If not found, append new element.
    // Refinement: Separate the search phase from the update/append action phase.
    // This highlights which parts are potentially parallelizable (search)
    // and which parts require synchronization (append, modifying *nzv).

    int k_found = -1; // Use a sentinel value to indicate if the index 'i' was found.
                      // Sparse array indices k are assumed 1-based and positive, so -1 is safe.

    // --- Data Dependency and Synchronization Consideration ---
    // The search loop reads from 'iv' and 'v' (implicitly, via k_found used later).
    // The append operation modifies '*nzv', 'v', and 'iv'. This is the critical section
    // if multiple threads update the *same* sparse vector concurrently.
    // The update operation (v[k_found] = val) is generally safe if index 'i' is unique
    // in 'iv' and multiple threads aren't trying to update the same k_found with different values.

    // Capture the current size of the sparse vector for the search loop bound.
    // This avoids potential issues if *nzv is modified by another thread concurrently
    // trying to append (in a scenario where multiple threads update the same vector).
    int current_nzv_limit = *nzv;

    // --- Search Phase: Identify if 'i' exists and where ---
    // This loop iterates through the existing non-zero elements.
    // The iterations are largely independent reads. In a parallel context
    // (e.g., #pragma omp for), finding the first match and setting k_found
    // would require careful handling (e.g., atomic update or critical section for k_found,
    // or processing all iterations and then checking k_found).
    // The original code iterated through all elements even if a match was found early;
    // we keep that behavior here for functional equivalence in structure.
    for (int k = 1; k <= current_nzv_limit; k++) {
        if (iv[k] == i) {
            k_found = k; // Record the index where 'i' was found.
            // No 'break' here, mirroring the original loop structure.
        }
    }

    // --- Action Phase: Update or Append based on Search Result ---
    // This section performs the actual modification to the sparse vector.
    // It acts based on the 'k_found' determined in the search phase.
    // The 'else' block (appending a new element) is the primary area
    // requiring synchronization (e.g., #pragma omp critical) if multiple threads
    // are calling vecset on the *same* sparse vector instance and new elements
    // might be added concurrently.

    if (k_found != -1) {
        // Load Imbalance / Distribution: The update is a constant-time operation once k_found is known.
        // Synchronization Overhead: Updating v[k_found] is generally low-overhead if k_found is unique per index 'i'.
        // Data Dependencies: Depends on k_found from the search phase.
        v[k_found] = val; // Index 'i' was found. Update the value at the recorded position.
    } else {
        // Load Imbalance / Distribution: Appending involves incrementing *nzv and assigning two values, also constant time.
        // Synchronization Overhead: Modifying *nzv and adding elements REQUIRES protection (e.g., critical section)
        // if multiple threads are appending to the same vector concurrently. This is the main bottleneck for parallel append.
        // Data Dependencies: Depends on *nzv's current value and increments it. Writes to *nzv, v, iv.
        *nzv = *nzv + 1;    // Increment the count of non-zero elements. CRITICAL SECTION CANDIDATE.
        v[*nzv] = val;      // Add the new value.
        iv[*nzv] = i;       // Add the new index.
    }
}