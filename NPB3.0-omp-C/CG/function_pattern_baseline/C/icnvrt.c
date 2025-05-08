static int icnvrt(double x, int ipwr2) {
#pragma omp declare simd
    return ((int)(ipwr2 * x));
}