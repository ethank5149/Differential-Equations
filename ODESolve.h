#pragma once
#include <stdio.h>
#include <math.h>

int forwardeuler(double (*function_pointer)(double, double), double x, double* y, int N, double dx);