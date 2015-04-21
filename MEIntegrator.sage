from scipy.integrate import ode
from numpy import array

class MEIntegrator(object):
    """
    Integrates master equation

    """
    def __init__(self, rho, H, L, dt):
        self.rho = rho
        self.H = H
        self.L = L
        self.dt = dt
        self.matshape = array(rho).shape
        self.InitODE()

    def Flat(self, A):
        """
        Returns one-dimensional numpy array representation of matrix A

        """
        arr = array(A)
        vec = arr.ravel()
        return vec

    def Reshape(self, vec):
        """
        Converts one-dimensional numpy array to square matrix

        """
        arr = vec.reshape(self.matshape)
        mat = matrix(arr)
        return mat

    def InitODE(self):
        def f(t, y):
            """
            Defines right-hand side of master equation ode

            Can be used with scipy.integrate.ode ode solver

            """
            return self.RHS_array(y)

        y = self.Flat(self.rho)
        self.ode = ode(f)
        self.ode.set_integrator('zvode')
        self.ode.set_initial_value(y)

    def RHS(self, rho):
        """
        Defines right-hand side of master equation ode

        """
        unitary_term = CDF(-I) * self.H.commutator(rho)

        matrix_space = rho.parent()
        zero_matrix = matrix_space.matrix()

        lindblad_term = zero_matrix
        for L in self.L:
            Lc = L.conjugate_transpose()
            lindblad_term += L * rho * Lc - 0.5 * rho.anticommutator(Lc * L)

        return unitary_term + lindblad_term

    def RHS_array(self, rho_flat):
        """
        Defines right-hand side of master equation ode

        Accepts density matrix in flat format

        """
        rho = self.Reshape(rho_flat)
        rhs = self.RHS(rho)
        rhs_flat = self.Flat(rhs)
        return rhs_flat

    def Step(self):
        """
        Computes time evolution of rho over time dt

        """
        t0 = self.ode.t
        t1 = t0 + self.dt
        self.ode.integrate(t1)
        self.rho = self.Reshape(self.ode.y)

    def Integrate(self, steps):
        """
        Do 'steps' integration steps

        Return list of values on each step

        """
        rho_list = [self.rho]
        for t in range(1, steps):
            self.Step()
            rho_list.append(self.rho)

        return rho_list
