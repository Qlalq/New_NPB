static int ilog2(int n) {
    // The original function computes ceil(log2(n)) for n > 1, and 0 for n == 1.
    // The loop `while (nn < n)` has a loop-carried dependency on `nn` and `lg`,
    // making internal parallelization difficult for a single value of n.
    // To improve structure for potential later OpenMP use (typically when
    // calling this function on many independent data points), we can
    // replace the data-dependent loop with a fixed sequence of steps
    // using bit manipulation. This approach is often faster and provides
    // a predictable control flow structure.

    // Handle the base case n = 1
    if (n == 1) {
        return 0;
    }

    // For n > 1, we want to compute ceil(log2(n)).
    // This is equivalent to floor(log2(n - 1)) + 1 for n > 1.
    // We will compute floor(log2(v)) for v = n - 1 using bitwise operations.
    // This method finds the position of the most significant bit in v,
    // which corresponds to floor(log2(v)) for v > 0.

    // Using unsigned int for bitwise operations is safer.
    unsigned int v = n - 1;

    // Initialize result for floor(log2(v)).
    int result_floor_log2 = 0;

    // The following stages perform a binary search on the bit position
    // to find the highest set bit. This replaces the iterative doubling
    // with a fixed sequence of checks and shifts.

    // Stage 1: Check if the MSB is in the upper 16 bits (for 32-bit int)
    if (v >= (1U << 16)) {
        result_floor_log2 += 16;
        v >>= 16; // Reduce problem to the upper 16 bits
    }

    // Stage 2: Check if the MSB is in the next 8 bits
    if (v >= (1U << 8)) {
        result_floor_log2 += 8;
        v >>= 8; // Reduce problem to the next 8 bits
    }

    // Stage 3: Check if the MSB is in the next 4 bits
    if (v >= (1U << 4)) {
        result_floor_log2 += 4;
        v >>= 4; // Reduce problem to the next 4 bits
    }

    // Stage 4: Check if the MSB is in the next 2 bits
    if (v >= (1U << 2)) {
        result_floor_log2 += 2;
        v >>= 2; // Reduce problem to the next 2 bits
    }

    // Stage 5: Check if the MSB is in the last 2 bits (values 2 or 3)
    if (v >= (1U << 1)) { // If v is 2 or 3
        result_floor_log2 += 1;
        // No need to shift further, v is now 0 or 1
    }

    // result_floor_log2 now holds floor(log2(n - 1)) for n > 1.
    // The original function returned ceil(log2(n)) for n > 1,
    // which is floor(log2(n - 1)) + 1 for n > 1.

    return result_floor_log2 + 1;
}