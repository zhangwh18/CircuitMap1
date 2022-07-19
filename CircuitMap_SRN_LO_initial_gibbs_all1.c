#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define Pi 3.1415926535897932384626433832795
#define M 1
#define g 1
#define T pow(10, 7)
#define N 5000
#define theta 0.5
double complex x0A[N], x0B[N];
// double complex x0AM[N], x0BM[N];
double complex phiA[N], phiB[N];
double Meanvalue[N], Sumvalue[N];
double Variance[1000];
int taxis[1000];
int main()
{

    double x, y, phi;

    // printf("%f", T);
    clock_t start_run, finish_run;
    double time_count;
    int n, m, t, Mind, taxisind = 0;
    double tsave = 1, sumMeanvalue, averMeanvalue, sumAvalue, sumBvalue;

    sumAvalue = 0;
    sumBvalue = 0;
    for (n = 0; n < N; n++)
    {
        y = (double)rand() / (double)((unsigned)RAND_MAX + 1);
        x = -log(1 - y);
        phi = (double)rand() / (double)((unsigned)RAND_MAX + 1) * (2 * Pi);
        x0A[n] = x * exp(I * phi);
        y = (double)rand() / (double)((unsigned)RAND_MAX + 1);
        x = -log(1 - y);
        phi = (double)rand() / (double)((unsigned)RAND_MAX + 1) * (2 * Pi);
        x0B[n] = x * exp(I * phi);
        sumAvalue = sumAvalue + cabs(x0A[n]) * cabs(x0A[n]) + cabs(x0B[n]) * cabs(x0B[n]);
        // sumBvalue = sumBvalue + cabs(x0B[n]) * cabs(x0B[n]);
    }

    for (n = 0; n < N; n++)
    {
        x0A[n] = x0A[n] / pow(sumAvalue, 0.5);
        x0B[n] = x0B[n] / pow(sumAvalue, 0.5);
    }
    // printf("%f\n", cabs(x0A[m][1]));

    sumAvalue = 0;
    sumBvalue = 0;
    for (n = 0; n < N; n++)
    {
        sumAvalue = sumAvalue + cabs(x0A[n]) * cabs(x0A[n]) + cabs(x0B[n]) * cabs(x0B[n]); // sumBvalue = sumBvalue + cabs(x0B[n]) * cabs(x0B[n]);
    }
    printf("%f\n", sumAvalue);
    // printf("%f\n", sumBvalue);

    for (n = 0; n < N; n++)
    {
        Sumvalue[n] = cabs(x0A[n]);
    }
    for (t = 1; t < T; t++)
    {
        start_run = clock();
        //#omp parallel for

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
            phiA[n] = 0 + I * 0;
            phiB[n] = 0 + I * 0;
            if (n == 0)
            {
                // printf("%f\n", cabs(x0AM[1]));
                phiA[n] = cos(theta) * cos(theta) * x0A[n] - cos(theta) * sin(theta) * x0B[N - 1] + sin(theta) * sin(theta) * x0A[n + 1] + cos(theta) * sin(theta) * x0B[n];
                phiB[n] = cos(theta) * cos(theta) * x0B[n] - cos(theta) * sin(theta) * x0A[n] + sin(theta) * sin(theta) * x0B[N - 1] + cos(theta) * sin(theta) * x0A[n + 1];
            }
            else if (n == N - 1)
            {
                phiA[n] = cos(theta) * cos(theta) * x0A[n] - cos(theta) * sin(theta) * x0B[n - 1] + sin(theta) * sin(theta) * x0A[0] + cos(theta) * sin(theta) * x0B[n];
                phiB[n] = cos(theta) * cos(theta) * x0B[n] - cos(theta) * sin(theta) * x0A[n] + sin(theta) * sin(theta) * x0B[n - 1] + cos(theta) * sin(theta) * x0A[0];
            }
            else
            {
                // printf("%f\n", cabs(x0AM[n]));
                phiA[n] = cos(theta) * cos(theta) * x0A[n] - cos(theta) * sin(theta) * x0B[n - 1] + sin(theta) * sin(theta) * x0A[n + 1] + cos(theta) * sin(theta) * x0B[n];
                phiB[n] = cos(theta) * cos(theta) * x0B[n] - cos(theta) * sin(theta) * x0A[n] + sin(theta) * sin(theta) * x0B[n - 1] + cos(theta) * sin(theta) * x0A[n + 1];
            }
            // x0A[Mind][n] = exp(I * g * cabs(phiA) * cabs(phiA)) * phiA;
            // x0B[Mind][n] = exp(I * g * cabs(phiB) * cabs(phiB)) * phiB;
        }
        for (n = 0; n < N; n++)
        {
            x0A[n] = exp(I * g * cabs(phiA[n]) * cabs(phiA[n])) * phiA[n];
            x0B[n] = exp(I * g * cabs(phiB[n]) * cabs(phiB[n])) * phiB[n];
        }
        // printf("%f\n", cabs(x0A[2]));
        // sumAvalue = 0;
        // sumBvalue = 0;
        // for (n = 0; n < N; n++)
        // {
        //     sumAvalue = sumAvalue + cabs(x0A[n]) * cabs(x0A[n]) + cabs(x0B[n]) * cabs(x0B[n]); // sumBvalue = sumBvalue + cabs(x0B[n]) * cabs(x0B[n]);
        // }
        // printf("%f\n", sumAvalue + sumBvalue);
        // printf("%f\n", sumBvalue);

        for (n = 0; n < N; n++)
        {

            Sumvalue[n] = Sumvalue[n] + cabs(x0A[n]);
        }
        // printf("%f\n", cabs(x0A[2]));
        // printf("%f\n", Sumvalue[2] / (t + 1));
        if (t == tsave)
        {
            // printf("%d\n", t);
            tsave = (int)(ceil(1.021 * tsave));
            // printf("%f\n", ceil(1.021 * tsave));
            taxis[taxisind] = t;

            sumMeanvalue = 0;
            for (n = 0; n < N; n++)
            {
                Meanvalue[n] = Sumvalue[n] / (t + 1);
                sumMeanvalue = sumMeanvalue + cabs(Meanvalue[n]);
            }
            // printf("%f\n", cabs(Meanvalue[2]));
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

        // printf("%f\n", time_count);

        // printf("\n");
    }

    FILE *fp1;
    fp1 = fopen("HFt05_try_ini", "w");
    perror("fopen");

    for (n = 0; n < 1000; n++)
    {
        // printf("%d\t", taxis[n]);
        // printf("%.16f\n", Variance[n]);
        fprintf(fp1, "%d\t", taxis[n]);
        fprintf(fp1, "%.16f\n", Variance[n]);
    }
    fclose(fp1);

    return 0;
}
