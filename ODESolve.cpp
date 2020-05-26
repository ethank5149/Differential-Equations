#include <stdio.h>
#include <math.h>
#include "ODESolve.h"

int forwardeuler(double (*function_pointer)(double, double), double x, double *y, int N, double dx)
{
    int i = 0;

    while (i < N - 1)
    {
        x += dx;
        y[i + 1] = y[i] + dx * (*function_pointer)(x, y[i]);
        i++;
    }
    return 0;
}