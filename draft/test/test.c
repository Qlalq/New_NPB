#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

#define N 512

void matrix_multiply(int A[N][N], int B[N][N], int C[N][N], int mode) {
    struct timeval start, end;
    gettimeofday(&start, NULL);

    switch (mode) {
        case 0: // 串行 -O0
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    C[i][j] = 0;
                    for (int k = 0; k < N; k++) {
                        C[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
            break;
        case 1: // 串行 -O3
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    C[i][j] = 0;
                    for (int k = 0; k < N; k++) {
                        C[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
            break;
        case 2: // OpenMP
            #pragma omp parallel for
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    C[i][j] = 0;
                    for (int k = 0; k < N; k++) {
                        C[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
            break;
        case 3: // -fopenmp-simd
            // #pragma omp simd
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    C[i][j] = 0;
                    for (int k = 0; k < N; k++) {
                        C[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
            break;
        case 4: // -ftree-parallelize-loops
            // #pragma omp parallel for collapse(2)
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    C[i][j] = 0;
                    for (int k = 0; k < N; k++) {
                        C[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
            break;
    }

    gettimeofday(&end, NULL);
    printf("%d: %lf 秒\n", mode, (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0);
}

int main(int argc, char *argv[]) {
    int A[N][N], B[N][N], C[N][N];
    
    // Initialize matrices
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = rand() % 10;
            B[i][j] = rand() % 10;
        }
    }

    if (argc < 2) {
        printf("Usage: %s <mode>\n", argv[0]);
        printf("0: 串行 -O0\n"); // gcc -O0 -o test test.c               0.723754 秒
        printf("1: 串行 -O3\n"); //gcc -O3 -o test test.c                0.216421 秒
        printf("2: OpenMP\n"); //gcc -fopenmp -o test test.c             0.062825 秒
        printf("3: -fopenmp-simd\n");//gcc -fopenmp-simd -o test test.c  0.970867 秒
        printf("4: -ftree-parallelize-loops\n");//gcc -ftree-parallelize-loops=2 -O3 -o test test.c           0.238057 秒
        return 1;
    }

    int mode = atoi(argv[1]);
    matrix_multiply(A, B, C, mode);

    return 0;
}