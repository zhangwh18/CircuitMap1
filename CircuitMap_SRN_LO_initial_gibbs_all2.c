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
#define g 1
#define T pow(10, 8)
#define N 5000
// #define theta 0.5
double complex x0A[N], x0B[N];
// double complex x0AM[N], x0BM[N];
double complex phiA[N], phiB[N];
double Meanvalue[N], Sumvalue[N];
double Variance[1000];
int taxis[1000];
int main()
{
    double theta = 0.05;
    double x, y, phi;

    // printf("%f", T);
    clock_t start_run, finish_run;
    double time_count;
    int n, m, t, Mind, taxisind = 0;
    double tsave = 1, sumMeanvalue, averMeanvalue, sumAvalue;

    sumAvalue = 0;

    for (n = 0; n < N; n++)
    {
        y = (double)rand() / (double)((unsigned)RAND_MAX + 1);
        x = -log(1 - y);
        phi = (double)rand() / (double)((unsigned)RAND_MAX) * (2 * Pi);
        x0A[n] = x * (cos(phi) + I * sin(phi));
        y = (double)rand() / (double)((unsigned)RAND_MAX + 1);
        x = -log(1 - y);
        phi = (double)rand() / (double)((unsigned)RAND_MAX) * (2 * Pi);
        x0B[n] = x * (cos(phi) + I * sin(phi));
        // printf("%f\n", cimag(x0A[n]));
        sumAvalue = sumAvalue + cabs(x0A[n]) * cabs(x0A[n]) + cabs(x0B[n]) * cabs(x0B[n]);
        // printf("%f\n", cimag(I * sin(phi)));
        // sumBvalue = sumBvalue + cabs(x0B[n]) * cabs(x0B[n]);
    }
    y = (double)rand() / (double)((unsigned)RAND_MAX + 1);
    x = -log(1 - y);
    // printf("%f\n", x);
    // printf("%f\n", cimag(x0A[1]));

    // x0A[N - 1] = x0A[0];
    // x0B[N - 1] = x0B[0];
    // // printf("%d\n", N - 1);
    // sumAvalue = sumAvalue + cabs(x0A[N - 1]) * cabs(x0A[N - 1]) + cabs(x0B[N - 1]) * cabs(x0B[N - 1]);

    for (n = 0; n < N; n++)
    {
        x0A[n] = x0A[n] * sqrt(N * 2) / pow(sumAvalue, 0.5);
        x0B[n] = x0B[n] * sqrt(N * 2) / pow(sumAvalue, 0.5);
        // printf("%f\n", cimag(x0A[n]));
    }
    // printf("%f\n", cabs(x0A[m][1]));

    // sumAvalue = 0;
    // sumBvalue = 0;
    // for (n = 0; n < N; n++)
    // {
    //     sumAvalue = sumAvalue + cabs(x0A[n]) * cabs(x0A[n]) + cabs(x0B[n]) * cabs(x0B[n]); // sumBvalue = sumBvalue + cabs(x0B[n]) * cabs(x0B[n]);
    // }
    // printf("%f\n", sumAvalue);
    for (n = 0; n < N; n++)
    {
        Sumvalue[n] = cabs(x0A[n]);
        // SumvalueB[n] = cabs(x0B[n]);
    }
    // printf("%f\n", Sumvalue[900]);
    // printf("%f\n", SumvalueB[900]);
    start_run = clock();

    for (t = 1; t < T; t++)
    {
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
                // printf("%f\n", cabs(x0AM[n]));
                phiA[n] = cos(theta) * cos(theta) * x0A[n] - cos(theta) * sin(theta) * x0B[n - 1] + sin(theta) * sin(theta) * x0A[n - 1] + cos(theta) * sin(theta) * x0B[n];
                phiB[n] = cos(theta) * cos(theta) * x0B[n] - cos(theta) * sin(theta) * x0A[n] + sin(theta) * sin(theta) * x0B[n + 1] + cos(theta) * sin(theta) * x0A[n + 1];
            }
            // x0A[Mind][n] = exp(I * g * cabs(phiA) * cabs(phiA)) * phiA;
            // x0B[Mind][n] = exp(I * g * cabs(phiB) * cabs(phiB)) * phiB;
        }
        for (n = 0; n < N; n++)
        {
            // printf("%f\n", cimag(x0B[n]));
            x0A[n] = (cos(g * cabs(phiA[n]) * cabs(phiA[n])) + I * sin(g * cabs(phiA[n]) * cabs(phiA[n]))) * phiA[n];
            x0B[n] = (cos(g * cabs(phiB[n]) * cabs(phiB[n])) + I * sin(g * cabs(phiB[n]) * cabs(phiB[n]))) * phiB[n];
            // printf("%f\n", creal(x0A[n]));
        }
        // printf("%f\n", cimag(x0A[1]));
        // sumAvalue = 0;
        // // sumBvalue = 0;
        // for (n = 0; n < N; n++)
        // {
        //     sumAvalue = sumAvalue + cabs(x0A[n]) * cabs(x0A[n]) + cabs(x0B[n]) * cabs(x0B[n]); // sumBvalue = sumBvalue + cabs(x0B[n]) * cabs(x0B[n]);
        // }
        // printf("%f\n", sumAvalue);
        for (n = 0; n < N; n++)
        {
            Sumvalue[n] = Sumvalue[n] + cabs(x0A[n]);
            // SumvalueB[n] = SumvalueB[n] + cabs(x0B[n]);
        }
        // printf("%f\n", cabs(x0A[1]));
        if (t == tsave)
        {
            // printf("%d\n", t);
            tsave = (int)(ceil(1.021 * tsave));
            taxis[taxisind] = t;

            sumMeanvalue = 0;
            for (n = 0; n < N; n++)
            {
                Meanvalue[n] = Sumvalue[n] / (t + 1);
                // MeanvalueB[n] = SumvalueB[n] / (t + 1);

                sumMeanvalue = sumMeanvalue + Meanvalue[n];
            }
            // FILE *fp2;
            // char str[20];
            // char dstr[20];
            // char Mstr[20];
            // strcpy(str, "TmeanAB_");
            // sprintf(dstr, "%.3lf", theta);
            // strcat(str, dstr);
            // strcat(str, "_");
            // sprintf(Mstr, "%d", t);
            // strcat(str, Mstr);
            // fp2 = fopen(str, "w");
            // for (n = 0; n < N; n++)
            // {
            //     fprintf(fp2, "%.16f\t", Meanvalue[n]);
            //     fprintf(fp2, "%.16f\n", MeanvalueB[n]);

            //     // fprintf(fp1, "%.16f\n", Variance[n]);
            // }
            // fclose(fp2);

            // printf("%f\n", cabs(Meanvalue[1]));
            averMeanvalue = sumMeanvalue / N;
            // printf("%f\n", cabs(averMeanvalue));

            sumMeanvalue = 0;
            for (n = 0; n < N; n++)
            {

                // printf("%f\n", (cabs(Meanvalue[n]) - averMeanvalue) * (cabs(Meanvalue[n]) - averMeanvalue));

                sumMeanvalue = sumMeanvalue + (Meanvalue[n] - averMeanvalue) * (Meanvalue[n] - averMeanvalue);
            }
            // printf("%f\n", cabs(sumMeanvalue));
            Variance[taxisind] = sumMeanvalue / N;
            taxisind = taxisind + 1;
        }

        // printf("\n");
    }
    finish_run = clock();
    time_count = (double)(finish_run - start_run) / CLOCKS_PER_SEC;

    // printf("%f\n", time_count);
    FILE *fp1;
    fp1 = fopen("Variance_g_1_N_5000_M_1_T_8_s_005", "w");
    // perror("fopen");

    // FILE *fp2;
    // fp2 = fopen("HFt1", "w");
    // perror("fopen");

    for (n = 0; n < 1000; n++)
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
