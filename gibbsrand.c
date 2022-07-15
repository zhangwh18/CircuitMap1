#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define Pi 3.1415926535897932384626433832795
#define M 5000
#define g 1
#define T pow(10, 3)
#define N 100000
#define theta 0
// double complex x0A[M][N], x0B[M][N];
// // double complex x0AM[N], x0BM[N];
// double complex phiA[N], phiB[N];
// double Meanvalue[M];
// double Variance[1000];
// int taxis[1000];
int main()
{

    double x[N], y[N], phi[N];

    // printf("%f", T);
    clock_t start_run, finish_run;
    double time_count;
    int n, m, t, Mind, taxisind = 0;

    for (n = 0; n < N; n++)
    {
        y[n] = (double)rand() / (double)((unsigned)RAND_MAX + 1);
        x[n] = -log(1.0 - y[n]);
        phi[n] = (double)rand() / (double)((unsigned)RAND_MAX + 1) * (2 * Pi);
        printf("%f %f %f\n", x[n], y[n], phi[n]);
    }
}