#include <stdio.h> // Assuming boolean is defined elsewhere or replace with int
#ifndef boolean
#define boolean int
#define FALSE 0
#define TRUE 1
#endif

static void sparse_refactored(
    double a[],		// Output: sparse matrix values (reused for intermediate)
    int colidx[],	// Output: sparse matrix column indices (reused for intermediate)
    int rowstr[],	// Output: sparse matrix row pointers (CSR format) (reused temporarily)
    int n,          // Dimension of vector x
    int arow[],		// Input: triplet row indices (1-based)
    int acol[],		// Input: triplet column indices (1-based)
    double aelt[],	// Input: triplet values
    int firstrow,   // First row of the input matrix block (1-based)
    int lastrow,    // Last row of the input matrix block (1-based)
    double x[],		// Input/Output: vector (used as temporary per row)
    boolean mark[],	// Input/Output: marker array (used as temporary per row)
    int nzloc[],	// Temporary: indices of non-zeros found in x for a row (used per row)
    int nnza)       // Number of input non-zero elements
{
    int nrows;
    int i, j, k; // Loop indices
    int nza_input_idx; // Index for iterating through input triplets (arow, acol, aelt)
    int current_nza_output; // Global counter for output non-zeros in the final phase
    double xi; // Temporary value from x

    // --- Phase 0: Calculate effective number of rows being processed ---
    nrows = lastrow - firstrow + 1;
    // This step is simple and independent.

    // --- Phase 1: Initialize temporary structures for counting ---
    // original: initializes rowstr (used for counts) and mark (part 1)
    // The original rowstr is overloaded. Use a temporary array specifically for counts.
    int row_counts[nrows + 2]; // Temporary storage for counts per row (1-based indexing for rows 1..nrows+1)

    // Initialize mark (for potential use in a first pass, though original code structure is slightly different)
    // This initialization loop is parallelizable.
    for (j = 1; j <= n; j++) {
        mark[j] = FALSE;
    }
    // Initialize temporary counts array to zero. This loop is parallelizable.
    for (j = 0; j <= nrows + 1; j++) {
         row_counts[j] = 0;
    }

    // --- Phase 2: Count non-zeros per row using input triplets ---
    // original: for (nza = 1; nza <= nnza; nza++) { ... rowstr[j]++; }
    // Dependency: rowstr[j] is modified by multiple iterations (histogramming).
    // This loop fills row_counts. This phase can be parallelized using atomics on row_counts
    // or by using thread-local counts and reducing them afterwards.
    for (nza_input_idx = 1; nza_input_idx <= nnza; nza_input_idx++) {
        // Map input row index to the processed block's index space (1-based relative to firstrow)
        // and then to the row_counts array index (shifted by +1 as per original code's rowstr usage).
        j = (arow[nza_input_idx] - firstrow + 1) + 1;
        row_counts[j]++; // This update needs atomic protection or equivalent in parallel.
    }
    // Original code set rowstr[n+1] = 0 here, likely for the initial counts structure.
    // This is handled implicitly by bounds or specific prefix sum logic now.

    // --- Phase 3: Compute prefix sum on row_counts for intermediate row pointers ---
    // original: rowstr[1] = 1; for (j = 2; ...) rowstr[j] = rowstr[j] + rowstr[j-1];
    // Dependency: Strong loop-carried dependency. This phase is inherently sequential or requires a specialized parallel scan algorithm.
    int intermediate_row_pointers[nrows + 2]; // Temporary storage for pointers into a/colidx for the intermediate matrix

    intermediate_row_pointers[0] = 0; // Using 0-based index 0 for convenience, although array access is 1-based
    intermediate_row_pointers[1] = 1; // The first element of the first row (1-based index)

    for (j = 2; j <= nrows + 1; j++) {
        intermediate_row_pointers[j] = intermediate_row_pointers[j-1] + row_counts[j];
    }
    // At this point, intermediate_row_pointers[j] (for j=1..nrows+1) gives the 1-based
    // index in a/colidx where elements of (input) row j-1+firstrow should start.

    // --- Phase 4: Zero out space in output arrays 'a' and 'colidx' for intermediate data ---
    // original: for(j = 0;...) for(k = rowstr[j];...) a[k] = 0.0;
    // This loop uses pointers calculated in Phase 3 (intermediate_row_pointers)
    // to zero the memory that will be used for the intermediate sparse matrix.
    // The outer loop (j) iterates over rows. This outer loop is parallelizable
    // because each iteration deals with a distinct segment of the a/colidx arrays,
    // defined by intermediate_row_pointers[j] and intermediate_row_pointers[j+1].
    for(j = 1; j <= nrows; j++) { // Iterate over rows 1..nrows in the processed block (maps to intermediate_row_pointers indices 1..nrows)
         // Inner loop sets elements for row j to zero.
         for(k = intermediate_row_pointers[j]; k < intermediate_row_pointers[j+1]; k++) {
	       a[k] = 0.0; // Using output 'a' as temporary storage for intermediate values
               colidx[k] = 0; // Using output 'colidx' as temporary storage for intermediate indices
         }
    }


    // --- Phase 5: Populate 'a' and 'colidx' with input triplets using intermediate pointers ---
    // original: for (nza = 1; ...) { j = arow[nza] - firstrow + 1; k = rowstr[j]; ... rowstr[j]++; }
    // Dependency: rowstr[j] is used as a mutable index and incremented within the loop over input triplets.
    // To parallelize this loop over nza_input_idx, we need a thread-safe way to get the next position
    // within each row's block in a/colidx. This can be done with atomics on per-row counters
    // or by pre-calculating each triplet's final position.
    // Using a temporary array for current fill positions initialized from intermediate_row_pointers.
    int current_fill_pos[nrows + 2];
    // Initialization loop is parallelizable.
    for(j=0; j<=nrows+1; j++) {
        current_fill_pos[j] = intermediate_row_pointers[j];
    }

    // The loop over input triplets (nnza).
    // If parallelized (e.g., across nza_input_idx), the access current_fill_pos[j]++
    // for a given j would be a race condition if multiple triplets map to the same row j.
    // This requires synchronization (atomics) on current_fill_pos array elements.
    for (nza_input_idx = 1; nza_input_idx <= nnza; nza_input_idx++) {
        j = arow[nza_input_idx] - firstrow + 1; // Row index in the processed block (1-based)
        k = current_fill_pos[j]; // Get current position for this row (this read is safe if writes are atomic)
        a[k] = aelt[nza_input_idx]; // Place value (using output 'a' as temporary intermediate storage)
        colidx[k] = acol[nza_input_idx]; // Place column index (using output 'colidx' as temporary intermediate storage)
        current_fill_pos[j]++; // Increment position for this row (this update needs atomic protection)
    }
    // After this phase, a and colidx contain the intermediate sparse matrix data
    // for the block of rows firstrow..lastrow, stored contiguously.

    // --- Phase 6: Fix output rowstr array based on intermediate pointers ---
    // original: for (j = nrows; ...) rowstr[j+1] = rowstr[j]; rowstr[1] = 1;
    // Dependency: This phase transforms the intermediate_row_pointers (which were conceptually
    // used to fill a/colidx) into the final output rowstr format by shifting.
    // This loop has a loop-carried dependency (rowstr[j+1] depends on rowstr[j]). It is sequential.
    // The original code reuses the 'rowstr' variable for the final output pointers.
    for (j = nrows; j >= 1; j--) {
        rowstr[j+1] = intermediate_row_pointers[j]; // Shift pointers
    }
    rowstr[1] = 1; // Set the first pointer (1-based index)

    // --- Phase 7: Initialize x and mark for the main calculation phase ---
    // original: for (i = 1; i <= n; i++) { x[i] = 0.0; mark[i] = FALSE; }
    // Dependency: Initialization loop is parallelizable.
    for (i = 1; i <= n; i++) {
        x[i] = 0.0;
        mark[i] = FALSE;
    }

    // --- Phase 8: Main calculation and build final output sparse matrix ---
    // original: The main loop iterates over rows (j) of the intermediate matrix,
    // performs a computation (SpMV-like x = x + A*1), collects results per row
    // (using x, mark, nzloc as temporary per-row storage), and appends results
    // to the final output a, colidx arrays, updating a global counter (nza)
    // and the output rowstr.
    // Dependency: The update of the global output non-zero counter (current_nza_output)
    // and the final output rowstr array within the outer loop body makes the
    // outer loop sequential in the original code's logic.
    // The inner parts (SpMV-like computation per row) can be structured for parallelism
    // *within* a thread processing a row, and x/mark/nzloc arrays can be thread-private.

    current_nza_output = 0; // Global counter for the final output matrix non-zeros (1-based)
    // int jajp1 = rowstr[1]; // Start index for the first row in the intermediate matrix (using the fixed rowstr)
                           // This variable isn't strictly necessary if using rowstr[j] directly

    // Outer loop iterates over rows (j) of the intermediate matrix (stored in a/colidx
    // with pointers now correctly in 'rowstr' from Phase 6).
    // This loop is sequential as the final output assembly logic relies on global state updates.
    for (j = 1; j <= nrows; j++) {
        int nzrow = 0; // Counter for non-zeros in the result of the current row's computation (per row)
        // Note: x, mark, nzloc are used here as temporary storage *per row*.
        // In a parallel version, these would be thread-private arrays for each thread
        // processing one or more rows. The maximum size of nzloc and the needed
        // part of x/mark is 'n'.

        // Inner loop: Perform SpMV-like calculation for row j of the intermediate matrix
        // using the intermediate data stored in a/colidx and pointed to by the fixed rowstr.
        // x[i] = x[i] + a[k] accumulates results for the current row in the temporary vector x.
        // Indices k from rowstr[j] to rowstr[j+1]-1 access the intermediate matrix.
        for (k = rowstr[j]; k < rowstr[j+1]; k++) { // Uses the 'fixed' rowstr (Phase 6)
            i = colidx[k]; // Column index from intermediate colidx
            x[i] = x[i] + a[k]; // Accumulate result in temporary vector x[i] (thread-private x assumed)

            // Mark columns that have a non-zero contribution from the current row's calculation
            // This logic uses the temporary mark array (thread-private assumed).
            if ( mark[i] == FALSE && x[i] != 0.0) { // Check if x[i] became non-zero for the first time for this row
                mark[i] = TRUE;
                nzrow = nzrow + 1;
                nzloc[nzrow] = i; // Store column index of touched element in temporary nzloc array (thread-private assumed)
            }
        }

        // Process the non-zeros found for the current row (stored via nzloc)
        // and append them to the final output arrays (a, colidx).
        // Dependency: This loop iterates through nzloc (which was populated sequentially per row).
        // Writing to the global a, colidx arrays and updating the global current_nza_output
        // and the output rowstr[j+1] is sequential globally in the original structure.
        for (k = 1; k <= nzrow; k++) {
            i = nzloc[k]; // Column index of the resulting non-zero for this row
            mark[i] = FALSE; // Reset mark[i] for the next row's processing (thread-private assumed)
            xi = x[i];       // Get the accumulated value from temporary x[i] (thread-private assumed)
            x[i] = 0.0;      // Reset x[i] for the next row's processing (confirms x is row-local temp storage)

            if (xi != 0.0) {
                // Append the resulting non-zero element to the global output matrix
                current_nza_output = current_nza_output + 1; // Increment global non-zero counter for output
                a[current_nza_output] = xi;       // Store value in final output array 'a'
                colidx[current_nza_output] = i; // Store column index in final output array 'colidx'
                // Dependency: These writes modify global shared arrays 'a' and 'colidx'
                // Dependency: The increment of 'current_nza_output' is a global shared counter update.
            }
        }

        // Update the row pointer for the *next* row (j+1) in the *final* output matrix
        // based on the global non-zero counter.
        rowstr[j+1] = current_nza_output + rowstr[1]; // Assuming rowstr[1] remains 1 (1-based index base)
        // Dependency: Modifies the final output rowstr array based on the global current_nza_output.

        // Note: The original code updated jajp1 = rowstr[j+1] here. This is just updating
        // the loop boundary variable for the *next* iteration of the outer loop.
        // It's not a dependency that forces sequentiality of the outer loop itself,
        // but rather a pattern of traversing the intermediate matrix using its rowstr.
        // The real sequential dependency is the global update of current_nza_output and rowstr[j+1].
    }

    // The final total number of non-zeros in the output matrix is current_nza_output.
    // This is available after the loop finishes.
    // The final output rowstr array is built sequentially within the loop.
}