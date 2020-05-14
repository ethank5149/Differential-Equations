#include <iostream>
#include <cmath>
#include "ODESolve.h"

constexpr auto PI = 3.14159265359;
constexpr auto N = 10000;

double f(double x,double y)
{
    return (cos(x)+0.0);
}

int main()
{
    double (*function_pointer)(double,double) = &f;
    double x0, y0, dx;

    x0 = 0.0;
    y0 = 0.0;
    dx = 0.5*PI/N;

    double y[N];
    y[0] = y0;

    forwardeuler(function_pointer,x0,y,N,dx);

    printf("%f\n", y[N-1]);
    
    return 0;
}