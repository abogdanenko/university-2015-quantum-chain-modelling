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
        self.integrator = ode(f)
        self.integrator.set_integrator('zvode')
        self.integrator.set_initial_value(y)

    def RHS(self, rho):
        """
        Defines right-hand side of master equation ode

        """
        unitary_term = CDF(-I) * self.H.commutator(rho)

        L = self.L
        Lc = L.conjugate_transpose()
        lindblad_term = L * rho * Lc - 0.5 * rho.anticommutator(Lc * L)

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
        self.integrator.integrate(self.dt)

    def Rho():
        """
        Return current rho

        """
        y = self.integrator.y
        rho = self.Reshape(y)
        return rho
