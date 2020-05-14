#include <stdio.h>
#include <math.h>
#include "ODESolve.h"


int forwardeuler(double (*function_pointer)(double,double), double x, double* y,int N, double dx)
{
    int i = 0;

    while (i < N-1) {
        x+=dx;
        y[i + 1] = y[i] + dx*(*function_pointer)(x, y[i]);
        i++;
    }
    return 0;
}

int backwardeuler(double (*function_pointer)(double, double), double x, double* y, int N, double dx)
{
    int i = 0;

    while (i < N - 1) {
        x += dx;

        /**/
        double bisection(double (*function_pointer)(double), double a, double b, double allowed_error, int max_iterations)
        {
            int itr = 0;
            double c, yc, yb, ya;

            ya = a - y[i] - dx * (*function_pointer)(x, a);
            yb = b - y[i] - dx * (*function_pointer)(x, b);

            if (ya * yb > 0.0)
            {
                printf("The bracket either isn't surrounding the root, or it's surrounding an even number of them\n");
                return 0;
            }

            if (a > b)
            {
                double temp = a;
                a = b;
                b = temp;
            }

            do
            {
                itr++;
                c = (a + b) / 2;

                yc = c - y[i] - dx * (*function_pointer)(x, c);
                yb = b - y[i] - dx * (*function_pointer)(x, b);

                if (yc * yb > 0.0)
                    b = c;
                else
                    a = c;

                if (fabs(b - a) < allowed_error)
                {
                    return c;
                }

            } while (itr < max_iterations);

            printf("\nMaximum Number of Iterations Reached Without Converging to Provided Tolerance\n");
            return c;
        }





        /**/
        /*Call this after updating x to use the next iterates value*/
        y[i + 1] = bisection(yn1-y[i] - dx * (*function_pointer)(x, yn1);
        i++;
    }
    return 0;
}

