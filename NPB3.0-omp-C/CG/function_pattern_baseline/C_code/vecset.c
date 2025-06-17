static void vecset(
    int n,
    double v[],    
    int iv[],    
    int *nzv,
    int i,
    double val)
{
    int k;
    boolean set;
    set = FALSE;

    // Parallelize the loop with OpenMP
    #pragma omp parallel for shared(set, iv, v, nzv, i, val) private(k) reduction(|:set)
    for (k = 1; k <= *nzv; k++) {
        if (iv[k] == i) {
            v[k] = val;
            set = TRUE;
        }
    }

    // Check and update *nzv if needed (critical section to ensure thread-safety)
    if (set == FALSE) {
        #pragma omp atomic
        (*nzv) = *nzv + 1;
        v[*nzv] = val;
        iv[*nzv] = i;
    }
}