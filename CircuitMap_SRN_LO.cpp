#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <omp.h>
using namespace std;
#define Pi 3.1415926535897932384626433832795
#define M 5000
#define g 1
#define T 1000
#define N 5000
#define theta 0.5
// std::complex<double> x0A[M][N], x0B[M][N];
double complex x0A[M][N], x0B[M][N];

int main()
{
    // std::complex<double> x0A[M][N], x0B[M][N];
    double Meanvalue[M];
    double Variance[10000];
    int taxis[10000];
    // printf("%f", T);
    clock_t start_run, finish_run;
    double time_count;
    int n, m, t, taxisind = 0;
    double tsave = 0, sumMeanvalue, averMeanvalue, sumAvalue[M], sumBvalue[M];
    for (m = 0; m < M; m++)
    {
        sumAvalue[m] = 0;
        sumBvalue[m] = 0;
        for (n = 0; n < N; n++)
        {
            x0A[m][n] = rand() + I * rand();
            x0B[m][n] = rand() + I * rand();
            sumAvalue[m] = sumAvalue[m] + cabs(x0A[m][n]);
            sumBvalue[m] = sumBvalue[m] + cabs(x0B[m][n]);
        }
    }
    for (m = 0; m < M; m++)
    {
        for (n = 0; n < N; n++)
        {
            x0A[m][n] = x0A[m][n] / sumAvalue[m];
            x0B[m][n] = x0B[m][n] / sumBvalue[m];
        }
        // printf("%f\n", cabs(x0A[m][1]));
    }
    sumAvalue[2] = 0;

    for (n = 0; n < N; n++)
    {
        sumAvalue[2] = sumAvalue[2] + cabs(x0A[2][n]);
    }
    printf("%f\n", sumAvalue[2]);

    for (n = 0; n < M; n++)
    {
        Meanvalue[n] = cabs(x0A[n][0]);
    }
    double complex phiA[N], phiB[N];
    double complex x0AM[N], x0BM[N];

    for (t = 0; t < T; t++)
    {
        start_run = clock();
        //#omp parallel for

#pragma omp parallel for private(phiA, phiB, x0AM, x0BM) shared(x0A, x0B)

        // printf("Hello from process: %d\n", omp_get_thread_num());

        // printf("%lu\n", sizeof(phiA) / sizeof(phiA[0]));
        // #pragma omp for

        for (int Mind = 0; Mind <= M - 1; Mind++)
        {
            // printf("Hi from process: %d\n", omp_get_thread_num());
            // printf("%d\n", Mind);
            for (int z = 0; z < N; z++)
            {
                // printf("%lu\n", sizeof(phiA) / sizeof(phiA[0]));
                // printf("%d\n", z);
                phiA[z] = 0 + I * 0;
                phiB[n] = 0 + I * 0;
                x0AM[z] = x0A[Mind][z];
                x0BM[z] = x0B[Mind][z];
            }
            // printf("%f\n", creal(x0A[0][0]));
            for (int nomp = 0; nomp < N; nomp++)
            {
                if (nomp == 0)
                {
                    // printf("%f\n", cabs(x0AM[1]));
                    phiA[nomp] = cos(theta) * cos(theta) * x0AM[nomp] - cos(theta) * sin(theta) * x0BM[N - 1] + sin(theta) * sin(theta) * x0AM[nomp + 1] + cos(theta) * sin(theta) * x0BM[nomp];
                    phiB[nomp] = cos(theta) * cos(theta) * x0BM[nomp] - cos(theta) * sin(theta) * x0AM[nomp] + sin(theta) * sin(theta) * x0BM[N - 1] + cos(theta) * sin(theta) * x0AM[nomp + 1];
                }
                if (nomp == N - 1)
                {
                    phiA[nomp] = cos(theta) * cos(theta) * x0AM[nomp] - cos(theta) * sin(theta) * x0BM[nomp - 1] + sin(theta) * sin(theta) * x0AM[0] + cos(theta) * sin(theta) * x0BM[nomp];
                    phiB[nomp] = cos(theta) * cos(theta) * x0BM[nomp] - cos(theta) * sin(theta) * x0AM[nomp] + sin(theta) * sin(theta) * x0BM[nomp - 1] + cos(theta) * sin(theta) * x0AM[0];
                }
                if (nomp<N - 1 & nomp> 0)
                {
                    // printf("%f\n", cabs(x0AM[n]));
                    phiA[nomp] = cos(theta) * cos(theta) * x0AM[nomp] - cos(theta) * sin(theta) * x0BM[nomp - 1] + sin(theta) * sin(theta) * x0AM[nomp + 1] + cos(theta) * sin(theta) * x0BM[nomp];
                    phiB[nomp] = cos(theta) * cos(theta) * x0BM[nomp] - cos(theta) * sin(theta) * x0AM[nomp] + sin(theta) * sin(theta) * x0BM[nomp - 1] + cos(theta) * sin(theta) * x0AM[nomp + 1];
                }
            }
            for (int nomp1 = 0; nomp1 < N; nomp1++)
            {
                x0A[Mind][nomp1] = (cos(g * cabs(phiA[nomp1]) * cabs(phiA[nomp1])) + I * sin(g * cabs(phiA[nomp1]) * cabs(phiA[nomp1]))) * (creal(phiA[nomp1]) + I * cimag(phiA[nomp1]));
                x0B[Mind][nomp1] = (cos(g * cabs(phiB[nomp1]) * cabs(phiB[nomp1])) + I * sin(g * cabs(phiB[nomp1]) * cabs(phiB[nomp1]))) * (creal(phiB[nomp1]) + I * cimag(phiB[nomp1]));
                // x0A[Mind][n] = exp(I * g * abs(phiA[n]) * abs(phiA[n])) * phiA[n];
                // x0B[Mind][n] = exp(I * g * abs(phiB[n]) * abs(phiB[n])) * phiB[n];
            }
        }
        // #pragma omp critical
        // printf("hello\n");

        for (n = 0; n < M; n++)
        {
            Meanvalue[n] = (cabs(Meanvalue[n]) * t + cabs(x0A[n][0])) / (t + 1);
        }

        // printf("hi");
        // printf("%d\n", abs(n));
        if (t == ceil(1.021 * (tsave + 1)))
        {
            // printf("%d\n", t);
            tsave = 1.021 * (tsave + 1);
            taxis[taxisind] = t;
            sumMeanvalue = 0;
            for (n = 0; n < M; n++)
            {
                sumMeanvalue = sumMeanvalue + cabs(Meanvalue[n]);
            }
            // printf("%f\n", cabs(Meanvalue[1]));
            averMeanvalue = sumMeanvalue / M;
            // printf("%f\n", cabs(averMeanvalue));

            sumMeanvalue = 0;
            for (n = 0; n < M; n++)
            {
                // printf("%f\n", (cabs(Meanvalue[n]) - averMeanvalue) * (cabs(Meanvalue[n]) - averMeanvalue));

                sumMeanvalue = sumMeanvalue + (cabs(Meanvalue[n]) - averMeanvalue) * (cabs(Meanvalue[n]) - averMeanvalue);
            }
            // printf("%f\n", cabs(sumMeanvalue));
            Variance[taxisind] = sumMeanvalue * 10000;
            taxisind = taxisind + 1;
            // printf("%d\n", taxisind);
        }
        finish_run = clock();
        time_count = (double)(finish_run - start_run) / CLOCKS_PER_SEC;
        // printf("\n");
        printf("%f\n", time_count);

        // printf("\n");
    }

    FILE *fp1;
    fp1 = fopen("HFt05_para", "w");
    // perror("fopen");

    // FILE *fp2;
    // fp2 = fopen("HFt1", "w");
    // perror("fopen");

    for (n = 0; n < 1000; n++)
    {
        fprintf(fp1, "%d\t", taxis[n]);
        fprintf(fp1, "%.8f\n", Variance[n]);
    }
    fclose(fp1);
    // fclose(fp2);

    // FILE *fp;
    // fp = fopen("samrat.txt", "w");
    // perror("fopen");
    // for (n = 0; n < 10; n++)
    // {
    //   fprintf(fp, "%d\n", 1);
    // }
    // fclose(fp);

    return 0;
}
