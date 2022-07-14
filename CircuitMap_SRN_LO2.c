#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <omp.h>
#define Pi 3.1415926535897932384626433832795
// #define M 10
#define g 1
#define T pow(10, 7)
#define N 5000
#define theta 0
double complex x0A[N], x0B[N];
double complex x0AM[N], x0BM[N];
double complex phiA, phiB;
double Meanvalue[N];
double Variance[10000];
int taxis[10000];
int main()
{
    // printf("%f", T);
    clock_t start_run, finish_run;
    double time_count;
    int n, m, t, Mind, taxisind = 0;
    double tsave = 0, sumMeanvalue, averMeanvalue, sumAvalue, sumBvalue;
    // for (m = 0; m < M; m++)
    // {
    sumAvalue = 0;
    sumBvalue = 0;
    for (n = 0; n < N; n++)
    {
        x0A[n] = rand() + I * rand();
        x0B[n] = rand() + I * rand();
        sumAvalue = sumAvalue + cabs(x0A[n]);
        sumBvalue = sumBvalue + cabs(x0B[n]);
    }
    // }
    // for (m = 0; m < M; m++)
    // {
    for (n = 0; n < N; n++)
    {
        x0A[n] = x0A[n] / sumAvalue;
        x0B[n] = x0B[n] / sumBvalue;
    }
    // printf("%f\n", cabs(x0A[m][1]));
    // }
    // sumAvalue[2] = 0;

    // for (n = 0; n < N; n++)
    // {
    //     sumAvalue[2] = sumAvalue[2] + cabs(x0A[2][n]);
    // }
    // printf("%f", sumAvalue[2]);

    for (n = 0; n < N; n++)
    {
        // x0A[n][0] = 1 + I * 0;
        Meanvalue[n] = cabs(x0A[n]);
    }
    for (t = 0; t < T; t++)
    {
        start_run = clock();
        //#omp parallel for
        // for (Mind = 0; Mind < M; Mind++)
        // {
        // for (n = 0; n < N; n++)
        // {
        //   phiA = 0 + I * 0;
        //   phiB = 0 + I * 0;
        //   x0AM[n] = x0A[Mind][n];
        //   x0BM[n] = x0B[Mind][n];
        // }
        // printf("%f\n", creal(x0A[0][0]));
        for (n = 0; n < N; n++)
        {
            phiA = 0 + I * 0;
            phiB = 0 + I * 0;
            if (n == 0)
            {
                // printf("%f\n", cabs(x0AM[1]));
                phiA = cos(theta) * cos(theta) * x0A[n] - cos(theta) * sin(theta) * x0B[N - 1] + sin(theta) * sin(theta) * x0A[n + 1] + cos(theta) * sin(theta) * x0B[n];
                phiB = cos(theta) * cos(theta) * x0B[n] - cos(theta) * sin(theta) * x0A[n] + sin(theta) * sin(theta) * x0B[N - 1] + cos(theta) * sin(theta) * x0A[n + 1];
            }
            else if (n == N - 1)
            {
                phiA = cos(theta) * cos(theta) * x0A[n] - cos(theta) * sin(theta) * x0B[n - 1] + sin(theta) * sin(theta) * x0A[0] + cos(theta) * sin(theta) * x0B[n];
                phiB = cos(theta) * cos(theta) * x0B[n] - cos(theta) * sin(theta) * x0A[n] + sin(theta) * sin(theta) * x0B[n - 1] + cos(theta) * sin(theta) * x0A[0];
            }
            else
            {
                // printf("%f\n", cabs(x0AM[n]));
                phiA = cos(theta) * cos(theta) * x0A[n] - cos(theta) * sin(theta) * x0B[n - 1] + sin(theta) * sin(theta) * x0A[n + 1] + cos(theta) * sin(theta) * x0B[n];
                phiB = cos(theta) * cos(theta) * x0B[n] - cos(theta) * sin(theta) * x0A[n] + sin(theta) * sin(theta) * x0B[n - 1] + cos(theta) * sin(theta) * x0A[n + 1];
            }
            x0A[n] = exp(I * g * cabs(phiA) * cabs(phiA)) * phiA;
            x0B[n] = exp(I * g * cabs(phiB) * cabs(phiB)) * phiB;
        }
        // for (n = 0; n < N; n++)
        // {
        //   x0A[Mind][n] = exp(I * g * cabs(phiA[n]) * cabs(phiA[n])) * phiA[n];
        //   x0B[Mind][n] = exp(I * g * cabs(phiB[n]) * cabs(phiB[n])) * phiB[n];
        // }
        // }
        for (n = 0; n < N; n++)
        {
            Meanvalue[n] = (cabs(Meanvalue[n]) * t + cabs(x0A[n])) / (t + 1);
        }
        // printf("%f\n", cabs(x0A[1][0]));
        if (t == ceil(1.021 * (tsave + 1)))
        {
            // printf("%d\n", t);
            tsave = 1.021 * (tsave + 1);
            taxis[taxisind] = t;
            sumMeanvalue = 0;
            for (n = 0; n < N; n++)
            {
                sumMeanvalue = sumMeanvalue + cabs(Meanvalue[n]);
            }
            // printf("%f\n", cabs(Meanvalue[1]));
            averMeanvalue = sumMeanvalue / N;
            // printf("%f\n", cabs(averMeanvalue));

            sumMeanvalue = 0;
            for (n = 0; n < N; n++)
            {
                // printf("%f\n", (cabs(Meanvalue[n]) - averMeanvalue) * (cabs(Meanvalue[n]) - averMeanvalue));

                sumMeanvalue = sumMeanvalue + (cabs(Meanvalue[n]) - averMeanvalue) * (cabs(Meanvalue[n]) - averMeanvalue);
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
    fp1 = fopen("HFt0_try2", "w");
    // perror("fopen");

    // FILE *fp2;
    // fp2 = fopen("HFt1", "w");
    // perror("fopen");

    for (n = 0; n < 10000; n++)
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
