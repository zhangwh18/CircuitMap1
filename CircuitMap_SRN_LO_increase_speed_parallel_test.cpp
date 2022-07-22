#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <string.h>
#include <omp.h>
using namespace std;
#define Pi 3.1415926535897932384626433832795
#define g 1
#define T pow(10, 5)
#define N 5000

int main(int argc, char **argv)
{
    double theta = 0.5;
    double x, y, phi, time_count, tsave = 1, averMeanvalue, sumAvalue = 0, sumMeanvalue;
    std::complex<double> x0A[N], x0B[N], phiA[N], phiB[N];
    double Meanvalue[N], Sumvalue[N];
    double Variance[1000];
    int taxis[1000];

    clock_t start_run, finish_run;
    int taxisind = -1;
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
        sumAvalue = sumAvalue + abs(x0A[n]) * abs(x0A[n]) + abs(x0B[n]) * abs(x0B[n]);
    }

    for (int n = 0; n < N; n++)
    {
        x0A[n] = x0A[n] * sqrt(N * 2) / pow(sumAvalue, 0.5);
        x0B[n] = x0B[n] * sqrt(N * 2) / pow(sumAvalue, 0.5);
        Sumvalue[n] = abs(x0A[n]);
    }
    int thread_id = 0, nloops;
    start_run = clock();
    for (int t = 1; t < 2000; t++)
    {
//#omp parallel for
#pragma omp parallel shared(phiA, phiB) private(thread_id, nloops)
        {
            nloops = 0;
            // int thread_id = 0;
#pragma omp for

            // printf("%d", thread_id);
            for (int n = 0; n < N; n++)
            {
                ++nloops;
                phiA[n] = 0 + I * 0;
                phiB[n] = 0 + I * 0;
                if (n == 0)
                {
                    phiA[n] = cos(theta) * cos(theta) * x0A[n] - cos(theta) * sin(theta) * x0B[N - 1] + sin(theta) * sin(theta) * x0A[N - 1] + cos(theta) * sin(theta) * x0B[n];
                    phiB[n] = cos(theta) * cos(theta) * x0B[n] - cos(theta) * sin(theta) * x0A[n] + sin(theta) * sin(theta) * x0B[n + 1] + cos(theta) * sin(theta) * x0A[n + 1];
                }
                else if (n == N - 1)
                {
                    phiA[n] = cos(theta) * cos(theta) * x0A[n] - cos(theta) * sin(theta) * x0B[n - 1] + sin(theta) * sin(theta) * x0A[n - 1] + cos(theta) * sin(theta) * x0B[n];
                    phiB[n] = cos(theta) * cos(theta) * x0B[n] - cos(theta) * sin(theta) * x0A[n] + sin(theta) * sin(theta) * x0B[0] + cos(theta) * sin(theta) * x0A[0];
                }
                else
                {
                    phiA[n] = cos(theta) * cos(theta) * x0A[n] - cos(theta) * sin(theta) * x0B[n - 1] + sin(theta) * sin(theta) * x0A[n - 1] + cos(theta) * sin(theta) * x0B[n];
                    phiB[n] = cos(theta) * cos(theta) * x0B[n] - cos(theta) * sin(theta) * x0A[n] + sin(theta) * sin(theta) * x0B[n + 1] + cos(theta) * sin(theta) * x0A[n + 1];
                }
            }
            thread_id = omp_get_thread_num();
            printf("Thread %d performed the %d thloop.\n", thread_id, nloops);
            // thread_id = omp_get_thread_num();
            // printf("Thread %d performed the loop.\n", thread_id);
        }

        for (int n = 0; n < N; n++)
        {
            x0A[n] = (cos(g * abs(phiA[n]) * abs(phiA[n])) + I * sin(g * abs(phiA[n]) * abs(phiA[n]))) * (real(phiA[n]) + I * imag(phiA[n]));
            x0B[n] = (cos(g * abs(phiB[n]) * abs(phiB[n])) + I * sin(g * abs(phiB[n]) * abs(phiB[n]))) * (real(phiB[n]) + I * imag(phiB[n]));
            Sumvalue[n] = Sumvalue[n] + abs(x0A[n]);
        }

        if (t == tsave)
        {
            tsave = (int)(ceil(1.021 * tsave));
            taxisind = taxisind + 1;
            taxis[taxisind] = t;
            averMeanvalue = 0;

            for (int n = 0; n < N; n++)
            {
                Meanvalue[n] = Sumvalue[n] / (t + 1);
                averMeanvalue = averMeanvalue + Meanvalue[n] / N;
            }

            sumMeanvalue = 0;

            for (int n = 0; n < N; n++)
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
    fp1 = fopen("Variance_g_1_N_5000_M_1_T_5_s_05", "w");
    for (int n = 0; n < 1000; n++)
    {
        fprintf(fp1, "%d\t", taxis[n]);
        fprintf(fp1, "%.16f\n", Variance[n]);
    }
    fclose(fp1);
    return 0;
}
