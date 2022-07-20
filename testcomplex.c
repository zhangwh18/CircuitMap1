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
#define g 10
#define T pow(10, 7)
#define N 5000

int main()
{
    double complex a, b, c;
    double phi = Pi / 4;
    a = 2.1 * (cos(phi) + I * sin(phi));
    b = 2 + I;
    c = a * b;
    printf("%f\n", cimag(c));
}