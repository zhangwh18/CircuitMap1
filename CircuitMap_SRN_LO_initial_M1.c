#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <omp.h>
#include <string.h>
#include <unistd.h>
#define Pi 3.1415926535897932384626433832795
// #define M 1
#define g 1
#define T pow(10, 8)
#define N 5000
// #define theta 0

int main(int argc, char *argv[])
{
    double theta = 0.5;
    double complex x0A[N], x0B[N];
    double complex phiA[N], phiB[N];
    double Meanvaluesave[1001], Meanvalue;
    int taxis[1001];
    int opt, M;
    opt = getopt(argc, argv, ": if: lrx");
    sscanf(argv[1], "%d", &M);

    // printf("%f", T);
    clock_t start_run, finish_run;
    double time_count;
    int n, t, taxisind = 0;
    double tsave = 0, sumAvalue = 0, sumBvalue = 0;
    start_run = clock();

    for (n = 0; n < N; n++)
    {
        x0A[n] = rand() + I * rand();
        x0B[n] = rand() + I * rand();
        sumAvalue = sumAvalue + cabs(x0A[n]);
        sumBvalue = sumBvalue + cabs(x0B[n]);
    }

    for (n = 0; n < N; n++)
    {
        x0A[n] = x0A[n] / sumAvalue;
        x0B[n] = x0B[n] / sumBvalue;
    }

    // for (n = 0; n < M; n++)
    // {
    //     x0A[n][n] = 1 + I * 0;
    //     Meanvalue[n] = cabs(x0A[n][0]);
    // }
    for (t = 0; t < T; t++)
    {

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
            // x0A[n] = exp(I * g * cabs(phiA) * cabs(phiA)) * phiA;
            // x0B[n] = exp(I * g * cabs(phiB) * cabs(phiB)) * phiB;
        }
        for (n = 0; n < N; n++)
        {
            x0A[n] = exp(I * g * cabs(phiA[n]) * cabs(phiA[n])) * phiA[n];
            x0B[n] = exp(I * g * cabs(phiB[n]) * cabs(phiB[n])) * phiB[n];
        }
        // for (n = 0; n < M; n++)
        // {
        Meanvalue = (cabs(Meanvalue) * t + cabs(x0A[0])) / (t + 1);
        // }
        // printf("%f\n", cabs(x0A[1][0]));
        if (t == ceil(1.021 * (tsave + 1)))
        {
            // printf("%d\n", t);
            tsave = 1.021 * (tsave + 1);
            taxis[taxisind] = t;
            Meanvaluesave[taxisind] = Meanvalue;
            taxisind = taxisind + 1;
        }

        // printf("\n");
    }

    FILE *fp1;
    char str[20];
    char dstr[20];
    char Mstr[20];
    strcpy(str, "Tmean_");
    sprintf(dstr, "%.3lf", theta);
    strcat(str, dstr);
    strcat(str, "_");
    sprintf(Mstr, "%d", M);
    strcat(str, Mstr);
    fp1 = fopen(str, "w");
    for (n = 0; n < 1001; n++)
    {
        fprintf(fp1, "%d\t", taxis[n]);
        fprintf(fp1, "%.16f\n", Meanvaluesave[n]);
    }
    fclose(fp1);
    finish_run = clock();
    time_count = (double)(finish_run - start_run) / CLOCKS_PER_SEC;

    printf("%f\n", time_count);
    return 0;
}
