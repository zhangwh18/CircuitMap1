#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#define Pi 3.1415926535897932384626433832795
#define g 1
#define T pow(10, 5)
#define N 5000

int main()
{
    double theta = 0.5;
    double x, y, phi, time_count, ct2 = cos(theta) * cos(theta), st2 = sin(theta) * sin(theta), sct = cos(theta) * sin(theta), tsave = 1, averMeanvalue, sumAvalue = 0, sumMeanvalue;
    double complex x0A[N], x0B[N], x0Are[N], x0Aim[N], x0Bre[N], x0Bim[N];
    double complex phiA[N], phiB[N], phiAre[N], phiBre[N], phiAim[N], phiBim[N];
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
        sumAvalue = sumAvalue + cabs(x0A[n]) * cabs(x0A[n]) + cabs(x0B[n]) * cabs(x0B[n]);
    }

    for (int n = 0; n < N; n++)
    {
        x0A[n] = x0A[n] * sqrt(N * 2) / pow(sumAvalue, 0.5);
        x0B[n] = x0B[n] * sqrt(N * 2) / pow(sumAvalue, 0.5);
        Sumvalue[n] = cabs(x0A[n]);
        x0Are[n] = creal(x0A[n]);
        x0Aim[n] = cimag(x0A[n]);
        x0Bre[n] = creal(x0B[n]);
        x0Bim[n] = cimag(x0B[n]);
    }

    start_run = clock();
    for (int t = 1; t < 10000; t++)
    {
        //#omp parallel for
        for (int n = 0; n < N; n++)
        {
            phiAre[n] = 0;
            phiBre[n] = 0;
            phiAim[n] = 0;
            phiBim[n] = 0;

            if (n == 0)
            {
                phiAre[n] = ct2 * x0Are[n] - sct * x0Bre[N - 1] + st2 * x0Are[N - 1] + sct * x0Bre[n];
                phiAim[n] = ct2 * x0Aim[n] - sct * x0Bim[N - 1] + st2 * x0Aim[N - 1] + sct * x0Bim[n];
                phiBre[n] = ct2 * x0Bre[n] - sct * x0Are[n] + st2 * x0Bre[n + 1] + sct * x0Are[n + 1];
                phiBim[n] = ct2 * x0Bim[n] - sct * x0Aim[n] + st2 * x0Bim[n + 1] + sct * x0Aim[n + 1];
            }
            else if (n == N - 1)
            {
                phiAre[n] = ct2 * x0Are[n] - sct * x0Bre[n - 1] + st2 * x0Are[n - 1] + sct * x0Bre[n];
                phiAim[n] = ct2 * x0Aim[n] - sct * x0Bim[n - 1] + st2 * x0Aim[n - 1] + sct * x0Bim[n];

                phiBre[n] = ct2 * x0Bre[n] - sct * x0Are[n] + st2 * x0Bre[0] + sct * x0Are[0];
                phiBim[n] = ct2 * x0Bim[n] - sct * x0Aim[n] + st2 * x0Bim[0] + sct * x0Aim[0];
            }
            else
            {
                phiAre[n] = ct2 * x0Are[n] - sct * x0Bre[n - 1] + st2 * x0Are[n - 1] + sct * x0Bre[n];
                phiAim[n] = ct2 * x0Aim[n] - sct * x0Bim[n - 1] + st2 * x0Aim[n - 1] + sct * x0Bim[n];

                phiBre[n] = ct2 * x0Bre[n] - sct * x0Are[n] + st2 * x0Bre[n + 1] + sct * x0Are[n + 1];
                phiBim[n] = ct2 * x0Bim[n] - sct * x0Aim[n] + st2 * x0Bim[n + 1] + sct * x0Aim[n + 1];
            }
        }
        for (int n = 0; n < N; n++)
        {
            // x0A[n] = (cos(g * cabs(phiA[n]) * cabs(phiA[n])) + I * sin(g * cabs(phiA[n]) * cabs(phiA[n]))) * phiA[n];
            // x0B[n] = (cos(g * cabs(phiB[n]) * cabs(phiB[n])) + I * sin(g * cabs(phiB[n]) * cabs(phiB[n]))) * phiB[n];
            x0Are[n] = cos(g * (phiAre[n] * phiAre[n]) + (phiAim[n] * phiAim[n])) * phiAre[n] - sin(g * (phiAre[n] * phiAre[n]) + (phiAim[n] * phiAim[n])) * phiAim[n];
            x0Aim[n] = sin(g * (phiAre[n] * phiAre[n]) + (phiAim[n] * phiAim[n])) * phiAre[n] + cos(g * (phiAre[n] * phiAre[n]) + (phiAim[n] * phiAim[n])) * phiAim[n];
            x0Bre[n] = cos(g * (phiBre[n] * phiBre[n]) + (phiBim[n] * phiBim[n])) * phiBre[n] - sin(g * (phiBre[n] * phiBre[n]) + (phiBim[n] * phiBim[n])) * phiBim[n];
            x0Bim[n] = sin(g * (phiBre[n] * phiBre[n]) + (phiBim[n] * phiBim[n])) * phiBre[n] + cos(g * (phiBre[n] * phiBre[n]) + (phiBim[n] * phiBim[n])) * phiBim[n];
            Sumvalue[n] = Sumvalue[n] + sqrt((phiAre[n] * phiAre[n]) + (phiAim[n] * phiAim[n]));
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
