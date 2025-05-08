#include <stdio.h> // For printf

// Note: T_MAX and timer_read(int) are assumed to be defined elsewhere
// and accessible in this scope.

// Structure to hold data for a single timer entry that needs printing.
// Using a struct allows us to collect all necessary pieces of information
// independently per timer before the final output phase.
typedef struct {
    int index;          // The index of the timer (0 to T_MAX-1)
    const char* name;   // Pointer to the static string name for the timer
    double value;       // The value read from the timer
} TimerEntry;

// Refined function to print timer values, structured for potential OpenMP parallelization.
// The core idea is to separate the independent data collection phase
// from the dependent output phase.
static void print_timers(void) {
    // Timer name strings - kept within the function as in the original.
    // Assuming T_MAX is intended to be the number of timers handled,
    // and this array has at least T_MAX elements, or we handle indices > array size.
    const char *tstrings[] = {
        "          total ",
        "          setup ",
        "            fft ",
        "         evolve ",
        "       checksum ",
        "         fftlow ",
        "        fftcopy "
    };
    // In a real case, T_MAX might be defined based on the size of tstrings,
    // e.g., #define T_MAX (sizeof(tstrings) / sizeof(tstrings[0]))

    // Array to temporarily store the data for non-zero timer entries found.
    // The size T_MAX is the maximum possible number of such entries.
    // This array serves as intermediate storage to decouple data gathering from printing.
    TimerEntry valid_entries[T_MAX];

    // Counter for the number of valid entries actually stored in valid_entries.
    // This needs careful handling (e.g., synchronization) if incremented in parallel.
    int valid_count = 0;

    // Phase 1: Collect timer data. This loop iterates through all potential timers.
    // The operations inside the loop (reading timer, accessing string, checking value)
    // are largely independent for each iteration 'i'.
    // This loop is a primary candidate for OpenMP parallelization (#pragma omp parallel for).
    // The critical part for parallelization is writing to the shared 'valid_entries'
    // array and incrementing 'valid_count'.
    for (int i = 0; i < T_MAX; i++) {
        double value = timer_read(i); // Read the timer value (assumed safe for concurrent reads)

        if (value != 0.0) {
            // If the timer value is non-zero, store its data.
            // In a parallel execution (#pragma omp parallel for), this entire block
            // needs to be protected by a critical section (#pragma omp critical)
            // to ensure atomic update of valid_entries and valid_count.
            // Alternatively, thread-local storage could be used followed by a merge,
            // which is more scalable for high thread counts or many non-zero timers.
            // For this refinement, structuring the sequential code to show the storage step suffices.

            // Check bounds before storing to prevent buffer overflow, just in case T_MAX
            // is defined larger than the capacity of valid_entries or the intent.
            if (valid_count < T_MAX) {
                valid_entries[valid_count].index = i;
                // Safely access the string name from the tstrings array.
                // Assumes T_MAX is <= the number of elements in tstrings, or handles the case.
                if (i < sizeof(tstrings) / sizeof(tstrings[0])) {
                     valid_entries[valid_count].name = tstrings[i];
                } else {
                     // Provide a fallback name if the index exceeds the tstrings array size.
                     // This handles potential mismatches between T_MAX and tstrings dimension.
                     valid_entries[valid_count].name = "Index Out of Bounds"; // Placeholder name
                     // A more sophisticated approach might format a string like "Timer %d" here.
                }
                valid_entries[valid_count].value = value;
                valid_count++; // Increment the counter. This needs atomic protection (#pragma omp atomic)
                               // or inclusion within the critical section if the loop is parallel.
            }
            // else { // valid_count == T_MAX
                 // Optionally handle the case where more than T_MAX non-zero timers are found
                 // if valid_entries array is of fixed T_MAX size and can't store more.
                 // For this example, we silently cap at T_MAX entries.
            // }
        }
    }

    // Phase 2: Print the collected data.
    // This loop iterates through the valid entries collected in Phase 1.
    // Printing to standard output (stdout) is a dependent operation; order matters.
    // Therefore, this phase must be executed sequentially, typically by a single thread
    // (e.g., the master thread after the parallel region finishes, or within a #pragma omp master block).
    for (int j = 0; j < valid_count; j++) {
        // Use the data stored in the valid_entries array.
        // This printf call is outside the parallel candidate loop (Phase 1).
        printf("timer %2d(%16s) :%10.6f\n",
               valid_entries[j].index,
               valid_entries[j].name,
               valid_entries[j].value);
    }

    // This structure effectively decouples the data gathering (parallel candidate)
    // from the ordered output (sequential).
}