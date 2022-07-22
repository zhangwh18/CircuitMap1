#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <string.h>

#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define Pi 3.1415926535897932384626433832795
#define M 1
#define g 10
#define T pow(10, 7)
#define N 5000
caA2 = cabs(phiA[n]) * cabs(phiA[n]);
x0A[n] = (caA2 + I * sin(caA2)) * phiA[n];
caA2 = cabs(phiB[n]) * cabs(phiB[n]);
x0B[n] = (cos(caA2) + I * sin(caA2)) * phiB[n];
int main()
{
    clock_t start_run, finish_run;
    float time_count, a;
    start_run = clock();
    for (int i = 5; i--;)
    {
        // a = 1;
        printf("%d\n", i);
    }
    finish_run = clock();
    time_count = (double)(finish_run - start_run) / CLOCKS_PER_SEC;
    printf("%f\n", time_count);
}