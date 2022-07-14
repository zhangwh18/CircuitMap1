#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <omp.h>
using namespace std;
#define Pi 3.1415926535897932384626433832795
#define M 4000
#define g 1
#define T pow(10, 6)
#define N 5000
#define theta 0.5
// std::complex<double> x0A[M][N], x0B[M][N];
// double complex x0AM[N], x0BM[N];
std::complex<double> x0A[M][N], x0B[M][N];

int main()
{
    // std::complex<double> x0A[M][N], x0B[M][N];
    double Meanvalue[M];
    double Variance[10000];
    int taxis[10000];
    // printf("%f", T);
    clock_t start_run, finish_run;
    double time_count;
    int k, t, taxisind = 0;
    double tsave = 0, sumMeanvalue, averMeanvalue, sumAvalue[M], sumBvalue[M];
    // for (m = 0; m < M; m++)
    // {
    //     sumAvalue[m] = 0;
    //     sumBvalue[m] = 0;
    //     for (n = 0; n < N; n++)
    //     {
    //         x0A[m][n] = rand() + I * rand();
    //         x0B[m][n] = rand() + I * rand();
    //         sumAvalue[m] = sumAvalue[m] + cabs(x0A[m][n]);
    //         sumBvalue[m] = sumBvalue[m] + cabs(x0B[m][n]);
    //     }
    // }
    // for (m = 0; m < M; m++)
    // {
    //     for (n = 0; n < N; n++)
    //     {
    //         x0A[m][n] = x0A[m][n] / sumAvalue[m];
    //         x0B[m][n] = x0B[m][n] / sumBvalue[m];
    //     }
    //     // printf("%f\n", cabs(x0A[m][1]));
    // }
    // sumAvalue[2] = 0;

    // for (n = 0; n < N; n++)
    // {
    //     sumAvalue[2] = sumAvalue[2] + cabs(x0A[2][n]);
    // }
    // printf("%f", sumAvalue[2]);

    for (k = 0; k < M; k++)
    {
        x0A[k][k] = 1 + I * 0;
        Meanvalue[k] = abs(x0A[k][0]);
    }
    std::complex<double> phiA, phiB;

    for (t = 0; t < T; t++)
    {
        start_run = clock();
#pragma omp parallel private(phiA, phiB) shared(x0A, x0B)
        {
#pragma omp for

            for (int Mind = 0; Mind < M; Mind++)
            {
                // printf("Hi from process: %d\n", omp_get_thread_num());
                // printf("%d\n", Mind);
                for (int n = 0; n < N; n++)
                {
                    phiA = 0 + I * 0;
                    phiB = 0 + I * 0;
                    if (n == 0)
                    {
                        // printf("%f\n", cabs(x0AM[1]));
                        phiA = cos(theta) * cos(theta) * x0A[Mind][n] - cos(theta) * sin(theta) * x0B[Mind][N - 1] + sin(theta) * sin(theta) * x0A[Mind][n + 1] + cos(theta) * sin(theta) * x0B[Mind][n];
                        phiB = cos(theta) * cos(theta) * x0B[Mind][n] - cos(theta) * sin(theta) * x0A[Mind][n] + sin(theta) * sin(theta) * x0B[Mind][N - 1] + cos(theta) * sin(theta) * x0A[Mind][n + 1];
                    }
                    else if (n == N - 1)
                    {
                        phiA = cos(theta) * cos(theta) * x0A[Mind][n] - cos(theta) * sin(theta) * x0B[Mind][n - 1] + sin(theta) * sin(theta) * x0A[Mind][0] + cos(theta) * sin(theta) * x0B[Mind][n];
                        phiB = cos(theta) * cos(theta) * x0B[Mind][n] - cos(theta) * sin(theta) * x0A[Mind][n] + sin(theta) * sin(theta) * x0B[Mind][n - 1] + cos(theta) * sin(theta) * x0A[Mind][0];
                    }
                    else
                    {
                        // printf("%f\n", cabs(x0AM[n]));
                        phiA = cos(theta) * cos(theta) * x0A[Mind][n] - cos(theta) * sin(theta) * x0B[Mind][n - 1] + sin(theta) * sin(theta) * x0A[Mind][n + 1] + cos(theta) * sin(theta) * x0B[Mind][n];
                        phiB = cos(theta) * cos(theta) * x0B[Mind][n] - cos(theta) * sin(theta) * x0A[Mind][n] + sin(theta) * sin(theta) * x0B[Mind][n - 1] + cos(theta) * sin(theta) * x0A[Mind][n + 1];
                    }
                    x0A[Mind][n] = (cos(g * abs(phiA) * abs(phiA)) + I * sin(g * abs(phiA) * abs(phiA))) * (real(phiA) + I * imag(phiA));
                    x0B[Mind][n] = (cos(g * abs(phiB) * abs(phiB)) + I * sin(g * abs(phiB) * abs(phiB))) * (real(phiB) + I * imag(phiB));
                }
                // for (n = 0; n < N; n++)
                // {
                //   x0A[Mind][n] = exp(I * g * cabs(phiA[n]) * cabs(phiA[n])) * phiA[n];
                //   x0B[Mind][n] = exp(I * g * cabs(phiB[n]) * cabs(phiB[n])) * phiB[n];
                // }
                printf("%d\n", Mind);
            }

            // #pragma omp critical
        }
        for (int n = 0; n < M; n++)
        {
            Meanvalue[n] = (abs(Meanvalue[n]) * t + abs(x0A[n][0])) / (t + 1);
        }

        // printf("%f\n", cabs(x0A[1][0]));
        if (t == ceil(1.021 * (tsave + 1)))
        {
            // printf("%d\n", t);
            tsave = 1.021 * (tsave + 1);
            taxis[taxisind] = t;
            sumMeanvalue = 0;
            for (int n = 0; n < M; n++)
            {
                sumMeanvalue = sumMeanvalue + abs(Meanvalue[n]);
            }
            // printf("%f\n", cabs(Meanvalue[1]));
            averMeanvalue = sumMeanvalue / M;
            // printf("%f\n", cabs(averMeanvalue));

            sumMeanvalue = 0;
            for (int n = 0; n < M; n++)
            {
                // printf("%f\n", (cabs(Meanvalue[n]) - averMeanvalue) * (cabs(Meanvalue[n]) - averMeanvalue));

                sumMeanvalue = sumMeanvalue + (abs(Meanvalue[n]) - averMeanvalue) * (abs(Meanvalue[n]) - averMeanvalue);
            }
            // printf("%f\n", cabs(sumMeanvalue));
            Variance[taxisind] = sumMeanvalue;
            taxisind = taxisind + 1;
        }
        finish_run = clock();
        time_count = (double)(finish_run - start_run) / CLOCKS_PER_SEC;

        printf("%f\n", time_count);

        // printf("\n");
    }

    FILE *fp1;
    fp1 = fopen("HFt05_try_ini", "w");
    // perror("fopen");

    // FILE *fp2;
    // fp2 = fopen("HFt1", "w");
    // perror("fopen");

    for (int n = 0; n < 10000; n++)
    {
        fprintf(fp1, "%d\t", taxis[n]);
        fprintf(fp1, "%.16f\n", Variance[n]);
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
