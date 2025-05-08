static int icnvrt(double x, int ipwr2) {
    // Refinement for OpenMP Readiness:
    // This function is inherently suitable for parallel execution without significant modifications
    // due to its simple nature and lack of internal dependencies or shared state.
    // 1. Data Dependencies: The computation (ipwr2 * x and cast) depends only on the function's
    //    input parameters (x, ipwr2). There are no internal loops or dependencies on prior states
    //    or shared memory locations. Each call is independent.
    // 2. Load Imbalance: The amount of computation (a single multiplication and cast) is constant
    //    for every call to this function. There is no internal load imbalance.
    // 3. Synchronization Overhead: The function operates only on local variables (input parameters
    //    passed by value and intermediate computation results). It does not access global or static
    //    variables (unless ipwr2 is global, but its usage here is as a parameter),
    //    and thus requires no internal synchronization mechanisms like locks or atomics.
    //
    // The structure below makes the intermediate calculation step explicit using a local variable,
    // which does not change its parallel characteristics but adds clarity to the computation flow.

    // Calculate the product as a double
    double temp_product = (double)ipwr2 * x;

    // Cast the result to an integer
    int final_result = (int)temp_product;

    return final_result;
}