import numpy as np
import matplotlib.pyplot as plt

# m * x'' + c * x' + k * x = 0
# x'' + 2 * b * x' + w0 ^ 2 * x = 0
# b = c / (2 * m)
# w0 = sqrt(k / m)


def exact(t, ics, params):
    f0, df0 = ics
    b, w0, w, A = params
    if b < w0:
        w1 = np.sqrt(w0**2-b**2)  # Underdamped:b<w0
        return np.exp(-b*t)*(f0*np.cos(w1*t)+(b*f0+df0)*np.sin(w1*t)/w1)
    elif b == w0:
        return np.exp(-b*t)*(f0+(b*f0+df0)*t)
    elif b > w0:
        w2 = np.sqrt(b**2-w0**2)  # Overdamped:b>w0
        return np.exp(-b*t)*(f0*np.cosh(w2*t)+(b*f0+df0)*np.sinh(w2*t)/w2)
    return 0


def dexact(t, ics, params):
    f0, df0 = ics
    b, w0, w, A = params
    f0 = f0 + 0j
    df0 = df0 + 0j
    b = b + 0j
    w0 = w0 + 0j
    C = np.sqrt(w0**2*(b**2-1))
    a = (2*C*np.exp(t*(b*w0+C)))**(-1)
    ans = a*(C*df0*(np.exp(2*C*t)+1)+C**2*f0*(np.exp(2*C*t)-1) -
             b*w0*(df0+b*w0*f0)*(np.exp(2*C*t)-1))
    return np.real(ans)


def F(t, w, A):
    return A*np.sin(w*t)


def solver(ics, params, T=10, dt=0.25, method='FD'):
    f0, df0 = ics
    b, w0, w, A = params

    Nt = int(T/dt)
    T = Nt*dt
    f = np.zeros(Nt+1)
    t = np.linspace(0, T, Nt+1)

    f[0] = f0
    f[1] = df0*dt+f[0]
    if method == 'FD':
        for n in range(0, Nt-1):
            f[n+2] = dt**2*F(n*dt, w, A)+2*(1-b*dt) * \
                f[n+1]-(1-2*b*dt+(w0*dt)**2)*f[n]
    if method == 'BD':
        for n in range(0, Nt-1):
            f[n+2] = (dt**2*F(n*dt, w, A)+2*(1+b*dt) *
                      f[n+1]-f[n])/(1+2*b*dt+(w0*dt)**2)
    if method == 'CN':
        for n in range(0, Nt-1):
            f[n+2] = (dt**2*F(n*dt, w, A)+(2-(w0*dt)**2) *
                      f[n+1]-(1-b*dt)*f[n])/(1+b*dt)
    return f, t


def diff(f, T, dt, method='FD'):
    Nt = int(T/dt)
    T = Nt*dt
    df = np.zeros(np.shape(f)[0])

    if method == 'FD':
        for n in range(0, Nt-2):
            df[n] = (f[n+1]-f[n])/dt
    if method == 'BD':
        for n in range(Nt, 2, -1):
            df[n] = (f[n]-f[n-1])/dt
    if method == 'CN':
        for n in range(1, Nt-2):
            df[n] = (f[n+1]-f[n-1])/(2*dt)
    return df


def example_phase_plot_comparison(ics, params, T, dt, save_fig=False):
    b, w0, w, A = params
    Nt = int(T/dt)
    T = Nt*dt
    te = np.linspace(0, T, Nt+1)

    f1, t1 = solver(ics=ics, params=params, T=T, dt=dt, method='FD')
    f2, t2 = solver(ics=ics, params=params, T=T, dt=dt, method='BD')
    f3, t3 = solver(ics=ics, params=params, T=T, dt=dt, method='CN')
    fe = exact(te, ics=ics, params=params)
    df1 = diff(f1, T, dt, method='FD')
    df2 = diff(f2, T, dt, method='BD')
    df3 = diff(f3, T, dt, method='CN')
    dfe = dexact(te, ics=ics, params=params)

    plt.plot(f1, df1, 'r-', label='Phase Plot - FD')
    plt.plot(f2, df2, 'b-', label='Phase Plot - BD')
    plt.plot(f3, df3, 'g-', label='Phase Plot - CN')
    plt.plot(fe, dfe, 'k-', label='Phase Plot - Exact')

    plt.title('Phase Plot Comparison')
    plt.ylabel('du/dt')
    plt.xlabel('u')
    plt.legend()
    plt.grid()

    if save_fig:
        plt.savefig(f"figures/harmonic_oscillator_phase_plots_b={b}_w0={w0}_dt={dt}.png", dpi=250)
        plt.show()
    else:
        plt.show()


def example_plot_comparison(ics, params, T, dt, save_fig=False):
    b, w0, w, A = params

    Nt = int(T/dt)
    T = Nt*dt
    te = np.linspace(0, T, Nt+1)

    f1, t1 = solver(ics=ics, params=params, T=T, dt=dt, method='FD')
    f2, t2 = solver(ics=ics, params=params, T=T, dt=dt, method='BD')
    f3, t3 = solver(ics=ics, params=params, T=T, dt=dt, method='CN')
    fe = exact(te, ics=ics, params=params)

    plt.subplot(2, 1, 1)
    plt.plot(t1, f1, 'r-', label='u(t) - FD')
    plt.plot(t2, f2, 'g-', label='u(t) - BD')
    plt.plot(t3, f3, 'b-', label='u(t) - CN')
    plt.title('Finite Difference Scheme Comparison')
    plt.ylabel('U(t)')
    plt.xlabel('t')
    plt.legend()
    plt.grid()

    plt.subplot(2, 1, 2)
    plt.plot(t1, f1-fe, 'r--', label='u(t) Error - FD')
    plt.plot(t2, f2-fe, 'g--', label='u(t) Error - BD')
    plt.plot(t3, f3-fe, 'b--', label='u(t) Error - CN')
    plt.title('Finite Difference Scheme Error')
    plt.ylabel('Error(t)')
    plt.xlabel('t')
    plt.legend()
    plt.grid()

    plt.tight_layout()

    if save_fig:
        plt.savefig(f"figures/harmonic_oscillator_plot_comparison_b={b}_w0={w0}_dt={dt}.png", dpi=250)
        plt.show()
    else:
        plt.show()


def example_plot(ics, params, T, dt, save_fig=False):
    b, w0, w, A = params

    Nt = int(T/dt)
    T = Nt*dt

    f1, t1 = solver(ics=ics, params=params, T=T, dt=dt, method='FD')
    f2, t2 = solver(ics=ics, params=params, T=T, dt=dt, method='BD')
    f3, t3 = solver(ics=ics, params=params, T=T, dt=dt, method='CN')

    plt.plot(t1, f1, 'r-', label='u(t) - FD')
    plt.plot(t2, f2, 'g-', label='u(t) - BD')
    plt.plot(t3, f3, 'b-', label='u(t) - CN')
    plt.title('Finite Difference Schemes')
    plt.ylabel('U(t)')
    plt.xlabel('t')
    plt.legend()
    plt.grid()

    if save_fig:
        plt.savefig(f"figures/harmonic_oscillator_plot_b={b}_w0={w0}_dt={dt}.png", dpi=250)
        plt.show()
    else:
        plt.show()


# Undamped: b = 0
# Underdamped: 0 <= b < 0
# Critically damped: b = 1
# Overdamped: b > 1

dt, T = 0.1, 50
u0, du0 = 1, -1
b, w0 = 0.5, 3.0

w, A = 1.0, 0.0

ics = (u0, du0)
params = (b, w0, w, A)

# example_phase_plot_comparison(ics, params, T, dt, save_fig=True)
example_plot_comparison(ics, params, T, dt, save_fig=True)
# example_plot(ics, params, T, dt, save_fig=True)
