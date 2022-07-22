// compute the sum of two arrays in parallel
#include <stdio.h>
#include <omp.h>
#include <time.h>

#define N 100000
int main(void)
{
    clock_t start_run, finish_run;
    double time_count;
    start_run = clock();

    float a[N], b[N], c[N];
    int i;

    /* Initialize arrays a and b */
    for (i = 0; i < N; i++)
    {
        a[i] = i * 2.0;
        b[i] = i * 3.0;
    }

    /* Compute values of array c = a+b in parallel. */
#pragma omp parallel shared(a, b, c) private(i)
    {
        int count = 0;
#pragma omp for
        for (i = 0; i < N; i++)
        {
            c[i] = a[i] + b[i];
            // printf("%f\n", c[10]);
            count++;
        }
        printf("process: % d iterated %d times\n", omp_get_thread_num(), count);
    }
    finish_run = clock();
    time_count = (double)(finish_run - start_run) / CLOCKS_PER_SEC;
    printf("%f\n", time_count);
}