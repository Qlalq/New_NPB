#include <omp.h>

static int icnvrt(double x, int ipwr2) {
    int result;
    
    #pragma omp parallel
    {
        #pragma omp single
        result = (int)(ipwr2 * x);
    }
    
    return result;
}