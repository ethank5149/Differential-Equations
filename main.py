from Runge_Kutta import RKMethods
import sympy as sp

import numpy as np
from sympy import symbols, Function, Eq
from sympy import Rational as R

n = symbols("n", integer=True)
y, t, h, alpha = symbols(r"y_n t_n h \alpha")
f = Function("f")

yn1 = sp.symbols("y_{n+1}")

ystarn1 = sp.symbols("y^*_{n+1}")
en1 = sp.symbols("e_{n+1}")
k = sp.IndexedBase("k")

methods = [RKMethods.euler,
           RKMethods.midpoint,
           RKMethods.rk2_heun,
           RKMethods.rk2_ralston,
           RKMethods.rk3_kutta,
           RKMethods.rk3_heun,
           RKMethods.rk3_ralston,
           RKMethods.ssprk3,
           RKMethods.rk4,
           RKMethods.rk4_38,
           RKMethods.rk1_2_3f,
           RKMethods.rk1_2_2f,
           RKMethods.rk23,
           RKMethods.rkf23_4,
           RKMethods.rkf23_3,
           RKMethods.rk3_2_4f,
           RKMethods.rkf34_f1,
           RKMethods.rkf34_f2,
           RKMethods.rkck32,
           RKMethods.rkf45,
           RKMethods.rkck43,
           RKMethods.rkf45_f1,
           RKMethods.rk4_5_6,
           RKMethods.dopri54fm,
           RKMethods.rk5_4_7f_eq3,
           RKMethods.rk5_4_7f_eq1,
           RKMethods.rk5_4_7f_eq2,
           RKMethods.rkck54,
           RKMethods.dopri54fc,
           RKMethods.dopri54fs,
           RKMethods.dopri54m,
           RKMethods.rk56,
           RKMethods.rkf56,
           RKMethods.dopri65s,
           RKMethods.dopri65c,
           RKMethods.dopri65m,
           RKMethods.rk6_7_10_verner,
           RKMethods.rkf67,
           RKMethods.rk7_8_13_verner,
           RKMethods.rkf78,
           RKMethods.rk8_7_13m,
           RKMethods.rk8_9_16
           ]


def display(rk):
    matrix, weights, nodes = rk.matrix, rk.weights, rk.nodes
    s = len(weights)

    k_vector = sp.Matrix([h * f(t + nodes[j] * h, y + sum([matrix[j][i] * k[i] for i in range(0, j)])) if j > 0 else h * f(t, y) for j in range(s)])
    # print("k_vector")
    print(sp.latex(k_vector))

    y_n_plus_1 = y + h * sum([weights[i] * k[i] for i in range(0, s)])
    # print("y_n_plus_1")
    print(sp.latex(y_n_plus_1))

    embedded =  isinstance(rk, RKMethods.EmbeddedExplicitRKMethod)
    if embedded:
        correctorweights = rk.correctorweights
        y_star_n_plus_1 = y + h * sum([correctorweights[i-1] * k[i] for i in range(1, s)])
        err_n_plus_1 = sp.simplify(y_n_plus_1 - y_star_n_plus_1)

        # print("y_star_n_plus_1")
        print(sp.latex(y_star_n_plus_1))

        # print("err_n_plus_1")
        print(sp.latex(err_n_plus_1))

    print("\n\n\n")


for rk in methods:
    display(rk)
