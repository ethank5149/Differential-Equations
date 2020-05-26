import numpy as np
import matplotlib.pyplot as plt


def solver(C, a, T, dt, method='FD'):
    """
    Solve u'=-a*u, u(0)=C, for t in (0,T] with steps of dt.

    Args:
        C (float): Initial condition.
        a (float): Function parameter.
        T (float): Simulation time.
        dt (float): Time step.
        method (str): Finite difference method.

    Returns:
        u: Finite difference solution.
        t: Finite difference time grid.
    """

    Nt = int(T/dt)               # number of time intervals
    T = Nt*dt                    # adjust T to fit time step dt
    u = np.zeros(Nt+1)           # array of u[n] values
    t = np.linspace(0, T, Nt+1)  # time mesh

    u[0] = C                     # assign initial condition
    if method == 'FD':
        for n in range(0, Nt):   # n=0,1,...,Nt-1
            u[n+1] = u[n]*(1 - a*dt)
    if method == 'BD':
        for n in range(0, Nt):   # n=0,1,...,Nt-1
            u[n+1] = u[n]/(1 + dt*a)
    if method == 'CN':
        for n in range(0, Nt):   # n=0,1,...,Nt-1
            u[n+1] = u[n]*(1 - 0.5*a*dt)/(1 + 0.5*dt*a)
    return u, t


def solver_save(C, a, T, dt, method='FD', filename='solution.dat'):
    """
    Solve u'=-a*u, u(0)=C, for t in (0,T] with steps of dt.
    Minimum use of memory. The solution is stored in a file
    (with name 'filename') for later plotting.

    Args:
        C (float): Initial condition parameter.
        a (float): Function parameter.
        T (float): Simulation time.
        dt (float): Time step.
        method (str): Finite difference method.

    Returns:
        None.
    """

    u, t = solver(C=C, a=a, T=T, dt=dt, method=method)
    sol = np.column_stack((u, t))
    np.savetxt(filename, sol)
    return None


def solver_read(filename='solution.dat'):
    """
    Read finite difference solution from solver_save().

    Args:
        filename: Solution data file, default is 'solution.dat'.

    Returns:
        u: Finite difference solution.
        t: Finite difference time grid.
    """
    data = np.loadtxt(filename)
    t = data[:, 1]
    u = data[:, 0]
    return u, t


def u_exact(t, C, a):
    """
    Gives exact solution of equation.

    Args:
        C (float): Initial condition parameter.
        a (float): Function parameter.
        t (float): Solution time.

    Returns:
        u: The exact answer answer at time t.
    """
    return C*np.exp(-a*t)


def example_plot_numerical_vs_exact(C, a, T, dt, method):
    """
    Compare the numerical vs exact solution in a plot.
    """

    u, t = solver(C=C, a=a, T=T, dt=dt, method=method)

    t_e = np.linspace(0, T, 1001)         # fine mesh for u_e
    u_e = u_exact(t_e, C, a)

    plt.plot(t, u, 'r--o',                # red dashes w/circles
             t_e, u_e, 'b-')              # blue line for exact sol.
    plt.legend(['Numerical', 'Exact'])
    plt.xlabel('t')
    plt.ylabel('u')
    plt.title(f'Method = {method}, dt = {dt}')
    plt.grid()
    plt.savefig(f'example_plot_numerical_vs_exact_{method}_{dt}.png')


def example_plot_method_comparison(C, a, T, dt):
    """
    Compare the different numerical solutions and exact solution all in one
    plot.
    """

    solver_save(C=C, a=a, T=T, dt=dt, method='FD', filename='solution_FD.dat')
    solver_save(C=C, a=a, T=T, dt=dt, method='BD', filename='solution_BD.dat')
    solver_save(C=C, a=a, T=T, dt=dt, method='CN', filename='solution_CN.dat')

    u_FD, t_FD = solver_read(filename='solution_FD.dat')
    u_BD, t_BD = solver_read(filename='solution_BD.dat')
    u_CN, t_CN = solver_read(filename='solution_CN.dat')

    t_e = np.linspace(0, T, 1001)        # fine mesh for u_e
    u_e = u_exact(t_e, C, a)

    plt.plot(
        t_FD,   u_FD,   'r--',            # red dashes w/circles
        t_BD,   u_BD,   'g--',            # green dashes w/circles
        t_CN,   u_CN,   'm--',            # magenta dashes w/circles
        t_e, u_e, 'b-')              # blue line for exact sol.
    plt.legend(['Numerical - FD', 'Numerical - BD', 'Numerical - CN', 'Exact'])
    plt.xlabel('t')
    plt.ylabel('u')
    plt.title(f'Numerical Method Comparison, dt={dt}')
    plt.grid()
    plt.savefig(f'example_plot_method_comparison_{dt}.png')


def example_output_method_error(C, a, T, dt, norm=2, methods=('FD', 'BD', 'CN')):
    """
    Run a case with the solver, compute error measure,
    and plot the numerical and exact solutions (if makeplot=True).
    """
    try:
        assert isinstance(norm, int)
    except AssertionError:
        print('Error: norm must be an integer')
        return None

    print('|method|norm|dt    |Error     |')
    for method in methods:
        u, t = solver(C, a, T, dt, method=method)    # Numerical solution
        u_e = u_exact(t, C, a)
        e = u_e - u
        E = (dt*sum(e**norm))**(1.0/norm)
        print(f'|{method:^6}|{norm:^4}|{dt:>1.4f}|{E:>1.4E}|')
    return None


# example_plot_method_comparison(C=10, a=2, T=15, dt=0.74)
# plt.show()


# example_plot_numerical_vs_exact(method='FD', C=1, a=2, T=10, dt=0.4)
# plt.show()


example_output_method_error(C=1, a=2, T=5, dt=0.25)
