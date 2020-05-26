# Author: Ethan Knox

'''
'''
import numpy as np
from matplotlib import pyplot as plt


class NDSolve:
    def __init__(self, f, U0, domain, max_iter=10, eps_iter=1e-7, rtol=1e-4, atol=1e-6, first_step='Default', min_step='Default', max_step='Default'):
        self.solver_methods = ['ForwardEuler', 'Leapfrog', 'Heun', 'RK2', 'RK4', 'RK3', 'AdamsBashforth2', 'AdamsBashforth3',
                               'AdamsBashMoulton2', 'AdamsBashforth4', 'AdamsBashMoulton3', 'EulerCromer', 'StaggeredMidpoint', 'MidpointIter', 'RKFehlberg']

        # Test first if U0 is sequence (len(U0) possible),
        # and use that as indicator for system of ODEs.
        try:
            self.neq = len(U0)
            U0 = np.asarray(U0)          # (assume U0 is sequence)
        except TypeError:
            # U0 has no __len__ method, assume it is a scalar
            self.neq = 1
            if isinstance(U0, int):
                U0 = float(U0)           # avoid integer division

        self.U0 = U0
        self.f = lambda u, t: np.asarray(f(u, t))
        self.t = np.asarray(domain)
        self.num_iterations = self.t.size - 1  # no of intervals

        if self.neq == 1:  # scalar ODEs
            self.u = np.zeros(self.num_iterations+1)
        else:              # systems of ODEs
            self.u = np.zeros((self.num_iterations+1, self.neq))
        self.u[0] = self.U0

        self.first_step = self.t[1] - \
            self.t[0] if first_step == 'Default' else first_step
        time_steps = self.t[1:] - self.t[:-1]
        self.min_step = 0.01*time_steps.min() if min_step == 'Default' else min_step
        self.max_step = 10*time_steps.max() if max_step == 'Default' else max_step
        self.rtol = rtol
        self.atol = atol
        self.max_iter = max_iter
        self.eps_iter = eps_iter

    def ForwardEuler(self):
        """
        Forward Euler scheme::

        u[n+1] = u[n] + dt*f(u[n], t[n])
        """
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            u, f, n, t = self.u, self.f, self.n, self.t
            dt = t[n+1] - t[n]
            u_new = u[n] + dt*f(u[n], t[n])
            self.u[n+1] = u_new
        return self.u, self.t

    def Leapfrog(self):
        """
        Leapfrog scheme::
            u[n+1] = u[n-1] + dt2*f(u[n], t[n])

        with::
            dt2 = t[n+1] - t[n-1]

        Forward Euler is used for the first step.
        """
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            u, f, n, t = self.u, self.f, self.n, self.t

            if n >= 1:
                dt2 = t[n+1] - t[n-1]
                u_new = u[n-1] + dt2*f(u[n], t[n])
            else:
                dt = t[n+1] - t[n]
                u_new = u[n] + dt*f(u[n], t[n])
            self.u[n+1] = u_new

        return self.u, self.t

    def Heun(self):
        """
        Heun's method, also known as an RungeKutta2 or Trapezoidal method.
        Basically, it is a central difference method, with one
        iteration and the Forward Euler scheme as start value.
        In this sense, it is a predictor-corrector method.

        Scheme::

            u[n+1] = u[n] + 0.5*dt* \
                (f(u[n],t[n]) + f(u[n]+dt*f(u[n],t[n]),t[n+1]))
        """
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            u, f, n, t = self.u, self.f, self.n, self.t
            dt = t[n+1] - t[n]
            u_star = u[n] + dt*f(u[n], t[n])  # Forward Euler step
            u_new = u[n] + 0.5*dt*(f(u[n], t[n]) + f(u_star, t[n+1]))
            self.u[n+1] = u_new

        return self.u, self.t

    def RK2(self):
        """
        Standard Runge-Kutta 2nd method::

            u[n+1] = u[n] + dt*f(u[n] + 0.5*(dt*f(u[n],t[n])),t[n] + 0.5*dt)
        """
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            u, f, n, t = self.u, self.f, self.n, self.t
            dt = t[n+1] - t[n]
            K1 = dt*f(u[n], t[n])
            K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
            u_new = u[n] + K2
            self.u[n+1] = u_new

        return self.u, self.t

    def RK4(self):
        """
        Standard RK4 method::

            u[n+1] = u[n] + (1/6.0)*(K1 + 2*K2 + 2*K3 + K4)

        where::
           K1 = dt*f(u[n], t[n])
           K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
           K3 = dt*f(u[n] + 0.5*K2, t[n] + 0.5*dt)
           K4 = dt*f(u[n] + K3, t[n] + dt)
        """
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            u, f, n, t = self.u, self.f, self.n, self.t
            dt = t[n+1] - t[n]
            dt2 = dt/2.0
            K1 = dt*f(u[n], t[n])
            K2 = dt*f(u[n] + 0.5*K1, t[n] + dt2)
            K3 = dt*f(u[n] + 0.5*K2, t[n] + dt2)
            K4 = dt*f(u[n] + K3, t[n] + dt)
            u_new = u[n] + (1/6.0)*(K1 + 2*K2 + 2*K3 + K4)
            self.u[n+1] = u_new

        return self.u, self.t

    def RK3(self):
        """
        RungeKutta3 method::

            u[n+1] = u[n] + (1/6.0)*(K1 + 4*K2 + K3)

        where::

            K1 = dt*f(u[n], t[n])
            K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
            K3 = dt*f(u[n] - K1 + 2*K2, t[n] + dt)
        """
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            u, f, n, t = self.u, self.f, self.n, self.t
            dt = t[n+1] - t[n]
            dt2 = dt/2.0
            K1 = dt*f(u[n], t[n])
            K2 = dt*f(u[n] + 0.5*K1, t[n] + dt2)
            K3 = dt*f(u[n] - K1 + 2*K2, t[n] + dt)
            u_new = u[n] + (1/6.0)*(K1 + 4*K2 + K3)
            self.u[n+1] = u_new

        return self.u, self.t

    def AdamsBashforth2(self):
        """
        Second-order Adams-Bashforth method::

            u[n+1] = u[n] + dt/2*(3*f(u[n], t[n]) - f(u[n-1], t[n-1]))

        for constant time step dt.

        RK2 is used for the first step.
        """
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            u, f, n, t = self.u, self.f, self.n, self.t

            if n >= 1:
                dt = t[n+1] - t[n]  # must be constant
                self.f_n = f(u[n], t[n])
                u_new = u[n] + dt/2.*(3*self.f_n - self.f_n_1)
                self.f_n_1 = self.f_n
            else:
                # Using RK2 for the first step
                dt = t[n+1] - t[n]
                K1 = dt*f(u[n], t[n])
                K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
                u_new = u[n] + K2
                self.f_n_1 = f(u[0], t[0])
            self.u[n+1] = u_new

        return self.u, self.t

    def AdamsBashforth3(self):
        """
        Third-order Adams-Bashforth method::

            u[n+1] = u[n] + dt/12.*(23*f(u[n], t[n]) - 16*f(u[n-1], t[n-1])
                                    + 5*f(u[n-2], t[n-2]))

        for constant time step dt.

        RK2 is used for the first steps.
        """
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            u, f, n, t = self.u, self.f, self.n, self.t

            if n >= 2:
                dt = t[n+1] - t[n]  # must be constant
                self.f_n = f(u[n], t[n])
                u_new = u[n] + dt/12. * \
                    (23*self.f_n - 16*self.f_n_1 + 5*self.f_n_2)
                self.f_n_1, self.f_n_2, self.f_n = self.f_n, self.f_n_1, self.f_n_2

            else:
                # Using RK2 for the first steps
                dt = t[n+1] - t[n]
                K1 = dt*f(u[n], t[n])
                K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
                u_new = u[n] + K2
                if n == 0:
                    self.f_n_2 = f(u[0], t[0])
                elif n == 1:
                    self.f_n_1 = f(u[1], t[1])
            self.u[n+1] = u_new

        return self.u, self.t

    def AdamsBashMoulton2(self):
        """
        Two-step (3rd-order) Adams-Bashforth method::

            predictor = u[n] + dt/12.*(23.*f(u[n], t[n]) - 16*f(u[n-1], t[n-1]) +
                                5*f(u[n-2], t[n-2]))
            corrector = u[n] + dt/12.*(8.*f(u[n], t[n]) - f(u[n-1], t[n-1]) +
                                5*f(predictor, t[n+1]))

        for constant time step dt.

        RK2 is used for the first two steps.
        """
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            u, f, n, t = self.u, self.f, self.n, self.t

            if n >= 2:
                dt = t[n+1] - t[n]  # must be constant
                self.f_n = f(u[n], t[n])
                predictor = u[n] + dt/12. * \
                    (23.*self.f_n - 16*self.f_n_1 + 5*self.f_n_2)
                u_new = u[n] + dt/12.*(8*self.f_n - self.f_n_1 +
                                       5*f(predictor, t[n + 1]))
                self.f_n_1, self.f_n_2 = self.f_n, self.f_n_1
            else:
                # Using RK2 for the first steps
                dt = t[n+1] - t[n]
                K1 = dt*f(u[n], t[n])
                K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
                u_new = u[n] + K2
                if n == 0:
                    self.f_n_2 = f(u[0], t[0])
                elif n == 1:
                    self.f_n_1 = f(u[1], t[1])

            self.u[n+1] = u_new

        return self.u, self.t

    def AdamsBashforth4(self):
        """
        Fourth-order Adams-Bashforth method::

            u[n+1] = u[n] + dt/24.*(55.*f(u[n], t[n]) - 59*f(u[n-1],
                                    t[n-1]) + 37*f(u[n-2], t[n-2]) - 9*f(u[n-3], t[n-3]))

        for constant time step dt.

        RK2 is used for the first three steps.
        """
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            u, f, n, t = self.u, self.f, self.n, self.t

            if n >= 3:
                dt = t[n+1] - t[n]  # must be constant
                self.f_n = f(u[n], t[n])
                u_new = u[n] + dt/24.*(55.*self.f_n - 59 *
                                       self.f_n_1 + 37*self.f_n_2 - 9*self.f_n_3)
                self.f_n_1, self.f_n_2, self.f_n_3 = self.f_n, self.f_n_1, self.f_n_2
            else:
                # Using RK2 for the first steps
                dt = t[n+1] - t[n]
                K1 = dt*f(u[n], t[n])
                K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
                u_new = u[n] + K2
                if n == 0:
                    self.f_n_3 = f(u[0], t[0])
                elif n == 1:
                    self.f_n_2 = f(u[1], t[1])
                elif n == 2:
                    self.f_n_1 = f(u[2], t[2])

            self.u[n+1] = u_new

        return self.u, self.t

    def AdamsBashMoulton3(self):
        """
        Three-step (4th-order) Adams-Bashforth method::

            predictor = u[n] + dt/24.*(55.*f(u[n], t[n]) - 59*f(
                u[n-1], t[n-1]) + 37*f(u[n-2], t[n-2]) - 9*f(u[n-3], t[n-3]))
            corrector = u[n] + dt/24.*(19.*f(u[n], t[n]) - 5*f(
                u[n-1], t[n-1]) + f(u[n-2], t[n-2]) + 9*f(predictor, t[n+1]))

        for constant time step dt.

        RK2 is used for the first three steps.
        """
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            u, f, n, t = self.u, self.f, self.n, self.t

            if n >= 3:
                dt = t[n+1] - t[n]  # must be constant
                self.f_n = f(u[n], t[n])
                predictor = u[n] + dt/24.*(55.*self.f_n - 59*self.f_n_1 +
                                           37*self.f_n_2 - 9*self.f_n_3)
                u_new = u[n] + dt/24.*(self.f_n_2 - 5*self.f_n_1 + 19*self.f_n +
                                       9*f(predictor, t[n + 1]))
                self.f_n_1, self.f_n_2, self.f_n_3 = self.f_n, self.f_n_1, self.f_n_2
            else:
                # Using RK2 for the first steps
                dt = t[n+1] - t[n]
                K1 = dt*f(u[n], t[n])
                K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
                u_new = u[n] + K2
                if n == 0:
                    self.f_n_3 = f(u[0], t[0])
                elif n == 1:
                    self.f_n_2 = f(u[1], t[1])
                elif n == 2:
                    self.f_n_1 = f(u[2], t[2])

            self.u[n+1] = u_new

        return self.u, self.t

    def EulerCromer(self):
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            u, f, n, t = self.u, self.f, self.n, self.t
            dt = t[n+1] - t[n]
            u_new = u[n].copy()
            # March forward velocities
            f_n = f(u[n], t[n])
            for i in range(0, len(u[n]), 2):
                u_new[i] = u[n, i] + dt*f_n[i]
            # March forward positions (u_new has now updated velocities, but old positions - but only the velocities are used in
            # the components of f(u_new, t) that enter the loop below, so the position values do not matter)
            f_np1 = f(u_new, t[n+1])
            for i in range(1, len(u[n]), 2):
                u_new[i] = u[n, i] + dt*f_np1[i]
            self.u[n+1] = u_new

        return self.u, self.t

    def StaggeredMidpoint(self):
        """
        Variant of the Euler-Cromer method based on staggered grid for
        positions and velocities and explicit centered (midpont) differences
        everywhere.
        Classical, intuitive, symplectic solver of second-order accuracy.

        u[0], u[2], etc. are velocities at t=(i+1/2)*dt.
        u[1], u[3], etc. are positions at t=i*dt.
        """
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            u, f, n, t = self.u, self.f, self.n, self.t
            dt = t[n+1] - t[n]
            u_new = u[n].copy()  # just need to allocate
            f_n = f(u[n], t[n])
            # March forward velocities
            if n == 0:
                # First step: use special formula for updating velocities
                for i in range(0, len(u[n]), 2):
                    u_new[i] = u[n, i] + dt/2*f_n[i]
            else:
                # Later: standard formula centered around n
                for i in range(0, len(u[n]), 2):
                    u_new[i] = u[n, i] + dt*f_n[i]
            # March forward positions (u_new has now updated velocities,
            # but old positions - but only the velocities are used in
            # the components of f(u_new, t) that enter the loop below,
            # so the position values do not matter)
            f_np1 = f(u_new, t[n+1])
            for i in range(1, len(u[n]), 2):
                u_new[i] = u[n, i] + dt*f_np1[i]
            self.u[n+1] = u_new

        return self.u, self.t

    def MidpointIter(self):
        """
        A midpoint/central difference method with max_iter fixed-point
        iterations to solve the nonlinear system.
        The Forward Euler scheme is recovered if max_iter=1 and f(u,t)
        is independent of t. For max_iter=2 we have the Heun/RK2 scheme.
        """

        # v is a help array needed in the method
        if self.neq == 1:
            # Scalar ODE: v can be one-dim array
            self.v = np.zeros(self.max_iter+1, self.u.dtype)
        else:
            # System of ODEs: v must be two-dim array
            self.v = np.zeros((self.max_iter+1, self.neq), self.u.dtype)
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            u, f, n, t, v = self.u, self.f, self.n, self.t, self.v
            dt = t[n+1] - t[n]

            v[0] = u[n]
            q = 0
            v_finished = False   # |v[q]-v[q-1]| < eps
            while not v_finished and q < self.max_iter:
                q += 1
                v[q] = u[n] + 0.5*dt*(f(v[q-1], t[n+1]) + f(u[n], t[n]))
                if abs(v[q] - v[q-1]).max() < self.eps_iter:
                    v_finished = True
                    self.num_iterations = q

            u_new = v[q]
            self.u[n+1] = u_new

        return self.u, self.t

    def RKFehlberg(self):
        """
        The classical adaptive Runge-Kutta-Fehlberg method of order 4-5.
        """

        # auxilatory function to pick up the middle number among 3 floats
        def middle(x, y=.1, z=4.):
            return sorted([x, y, z])[1]
        # The time loop
        self.n = 0  # time step counter
        for n in range(self.num_iterations):
            self.n = n
            f, n, rtol, atol = self.f, self.n, self.rtol, self.atol
            h, min_step, max_step = self.first_step, self.min_step, self.max_step
            u_n, t_n, t_np1 = self.u[n], self.t[n], self.t[n+1]

            dt = t_np1 - t_n

            # coefficients in Butcher tableau
            c = (1/4.,
                 3/8.,
                 3/32.,
                 9/32.,
                 12/13.,
                 1932/2197.,
                 -7200/2197.,
                 7296/2197.,
                 439/216.,
                 -8.,
                 3680/513.,
                 -845/4104.,
                 1/2.,
                 -8/27.,
                 2.,
                 -3544/2565.,
                 1859/4104.,
                 -11/40.,
                 1/360.,
                 -128/4275.,
                 -2197/75240.,
                 1/50.,
                 2/55.,
                 25/216.,
                 1408/2565.,
                 2197/4104.,
                 -1/5.)

            # u_i and t_i hold intermediate steps between t_n and t_np1
            u_i = [u_n]
            t_i = [t_n]
            t = t_n

            while abs(t - t_n) < abs(t_np1 - t_n):
                u, t = u_i[-1], t_i[-1]

                # internal steps
                k1 = h*f(u, t)
                k2 = h*f(u+k1*c[0], t+h*c[0])
                k3 = h*f(u+k1*c[2]+k2*c[3], t+h*c[1])
                k4 = h*f(u+k1*c[5]+k2*c[6]+k3*c[7], t+h*c[4])
                k5 = h*f(u+k1*c[8]+k2*c[9]+k3*c[10]+k4*c[11], t+h)
                k6 = h*f(u+k1*c[13]+k2*c[14]+k3*c[15]+k4*c[16]+k5*c[17],
                         t+h*c[12])
                u_new = u + k1*c[23] + k3*c[24] + k4*c[25] + k5*c[26]

                # local error between 2 levels
                error = np.abs(k1*c[18] + k3*c[19] + k4*c[20] +
                               k5*c[21] + k6*c[22])
                tol = rtol*np.abs(u_new) + atol
                # Error factor = local-error/error-tolerance
                rms = error/tol
                rms_norm = np.sqrt((np.sum(rms*rms))/self.neq)

                # Close enough or step size can not be reduced any more
                if rms_norm <= 1. or h <= min_step:
                    u_i.append(u_new)
                    t_i.append(t+h)

                # prevent the error of dividing absolute zero
                error = np.asarray([(1e-16 if x == 0. else x) for x in error]) \
                    if self.neq > 1 else (1e-16 if error == 0. else error)

                # Factor to adjust next step size
                s = (tol/(2*error))**0.25
                # Factor should be in a reasonable range[0.1,4.0]
                s = min(map(middle, s)) if self.neq > 1 else middle(s)

                # Step size should be in range [min_step, max_step]
                h = middle(h*s, self.min_step, self.max_step)

                # h should be set to 't_np1-t_i[-1]' at the last intern step.
                h = min(h, t_np1-t_i[-1])

            self.u[n+1] = u_new

        return self.u, self.t
