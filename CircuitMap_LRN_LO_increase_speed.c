#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#define Pi 3.1415926535897932384
// #define g 0.5
// #define T pow(10, 7)
// #define theta 1.0367
#define N 5000

int main()
{
    float g = 0.5, time_count;
    double x, y, phi, caA2, averMeanvalue, sumAvalue = 0, sumMeanvalue;
    double complex x0A[N], x0B[N];
    double complex phiA[N], phiB[N];
    double Meanvalue[N], Sumvalue[N], Variance[1000];
    int taxis[1000];
    clock_t start_run, finish_run;
    int taxisind = -1, tsave = 1;
    for (int n = 0; n < N; n++)
    {
        y = (double)rand() / (double)((unsigned)RAND_MAX + 1);
        x = -log(1 - y);
        phi = (double)rand() / (double)((unsigned)RAND_MAX) * (2 * Pi);
        x0A[n] = x * (cos(phi) + I * sin(phi));
        y = (double)rand() / (double)((unsigned)RAND_MAX + 1);
        x = -log(1 - y);
        phi = (double)rand() / (double)((unsigned)RAND_MAX) * (2 * Pi);
        x0B[n] = x * (cos(phi) + I * sin(phi));
        sumAvalue = sumAvalue + cabs(x0A[n]) * cabs(x0A[n]) + cabs(x0B[n]) * cabs(x0B[n]);
    }

    for (int n = 0; n < N; n++)
    {
        x0A[n] = x0A[n] * sqrt(N * 2) / pow(sumAvalue, 0.5);
        x0B[n] = x0B[n] * sqrt(N * 2) / pow(sumAvalue, 0.5);
        Sumvalue[n] = cabs(x0A[n]);
    }

    start_run = clock();
    for (int t = 1; t < 10000000; t++)
    {
        //#omp parallel for
        for (unsigned int n = N; n--;)
        {
            // phiA[n] = 0 + I * 0;
            // phiB[n] = 0 + I * 0;
            if (n == 0)
            {
                phiA[n] = 0.2591 * x0A[n] - 0.4382 * x0B[N - 1] + 0.7409 * x0A[N - 1] + 0.4382 * x0B[n];
                phiB[n] = 0.2591 * x0B[n] - 0.4382 * x0A[n] + 0.7409 * x0B[n + 1] + 0.4382 * x0A[n + 1];
            }
            else if (n == N - 1)
            {
                phiA[n] = 0.2591 * x0A[n] - 0.4382 * x0B[n - 1] + 0.7409 * x0A[n - 1] + 0.4382 * x0B[n];
                phiB[n] = 0.2591 * x0B[n] - 0.4382 * x0A[n] + 0.7409 * x0B[0] + 0.4382 * x0A[0];
            }
            else
            {
                phiA[n] = 0.2591 * x0A[n] - 0.4382 * x0B[n - 1] + 0.7409 * x0A[n - 1] + 0.4382 * x0B[n];
                phiB[n] = 0.2591 * x0B[n] - 0.4382 * x0A[n] + 0.7409 * x0B[n + 1] + 0.4382 * x0A[n + 1];
            }
        }
        for (unsigned int n = N; n--;)
        {
            caA2 = creal(phiA[n]) * creal(phiA[n]) + cimag(phiA[n]) * cimag(phiA[n]);
            // caA2 = caA2 * caA2;
            x0A[n] = (cos(0.5 * caA2) + I * sin(0.5 * caA2)) * phiA[n];
            caA2 = creal(phiB[n]) * creal(phiB[n]) + cimag(phiB[n]) * cimag(phiB[n]);
            // caA2 = caA2 * caA2;
            x0B[n] = (cos(0.5 * caA2) + I * sin(0.5 * caA2)) * phiB[n];
            // x0A[n] = (cos(cabs(phiA[n]) * cabs(phiA[n])) + I * sin(cabs(phiA[n]) * cabs(phiA[n]))) * phiA[n];
            // x0B[n] = (cos(cabs(phiB[n]) * cabs(phiB[n])) + I * sin(cabs(phiB[n]) * cabs(phiB[n]))) * phiB[n];
            // x0A[n] = cos(g * cabs(phiA[n]) * cabs(phiA[n])) * creal(phiA[n]) - sin(g * cabs(phiA[n]) * cabs(phiA[n])) * cimag(phiA[n]) + I * (sin(g * cabs(phiA[n]) * cabs(phiA[n])) * creal(phiA[n]) + cos(g * cabs(phiA[n]) * cabs(phiA[n])) * cimag(phiA[n]));
            // x0B[n] = cos(g * cabs(phiB[n]) * cabs(phiB[n])) * creal(phiB[n]) - sin(g * cabs(phiB[n]) * cabs(phiB[n])) * cimag(phiB[n]) + I * (sin(g * cabs(phiB[n]) * cabs(phiB[n])) * creal(phiB[n]) + cos(g * cabs(phiB[n]) * cabs(phiB[n])) * cimag(phiB[n]));
            // x0A[n] = cos(g * cabs(phiA[n]) ) * phiA[n] - sin(g * cabs(phiA[n])) * phiA[n] + I * (sin(g * cabs(phiA[n]) * *2) * phiA[n] + cos(g * cabs(phiA[n]) * *2) * phiA[n]);
            // x0B[n] = cos(g * cabs(phiB[n])) * phiB[n] - sin(g * (phiB[n] * phiB[n]) + (phiB[n] * phiB[n])) * phiB[n] + I * (sin(g * (phiB[n] * phiB[n]) + (phiB[n] * phiB[n])) * phiB[n] + cos(g * (phiB[n] * phiB[n]) + (phiB[n] * phiB[n])) * phiB[n]);
            Sumvalue[n] = Sumvalue[n] + cabs(x0A[n]);
        }
        if (t == tsave)
        {
            tsave = (int)(ceil(1.021 * tsave));
            taxisind = taxisind + 1;
            taxis[taxisind] = t;
            averMeanvalue = 0;
            for (unsigned int n = N; n--;)
            {
                Meanvalue[n] = Sumvalue[n] / (t + 1);
                averMeanvalue = averMeanvalue + Meanvalue[n] / N;
            }
            sumMeanvalue = 0;
            for (unsigned int n = N; n--;)
            {
                sumMeanvalue = sumMeanvalue + (Meanvalue[n] - averMeanvalue) * (Meanvalue[n] - averMeanvalue);
            }
            Variance[taxisind] = sumMeanvalue / N;
        }
    }
    finish_run = clock();
    time_count = (double)(finish_run - start_run) / CLOCKS_PER_SEC;
    printf("%f\n", time_count);
    FILE *fp1;
    fp1 = fopen("Variance_LRN_g_1_N_5000_M_1_T_7_g_05", "w");
    for (int n = 0; n < 1000; n++)
    {
        fprintf(fp1, "%d\t", taxis[n]);
        fprintf(fp1, "%.16f\n", Variance[n]);
    }
    fclose(fp1);
    return 0;
}
