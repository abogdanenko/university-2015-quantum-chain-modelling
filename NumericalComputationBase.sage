from scipy.integrate import ode
from numpy import array

class NumericalComputationBase(object):
    """
    Computes time evolution numerically

    """
    def __init__(self, sym):
        self.sym = sym
        self.params = NumericalParams()
        self.InitOperators()

    def IterationTime(self, t):
        """
        Returns moment of time corresponding to iteration number t

        """
        return t / self.params.time_steps * self.params.time_end

    def Subs(self, expr):
        """
        Returns expr with variables substituted with numbers

        """
        result = expr.subs(
            alpha = self.params.alpha,
            beta = self.params.beta,
            omega_a = self.params.omega_a,
            omega_c = self.params.omega_c)
        return result

    def SubsNum(self, expr):
        """
        Returns expr coerced to CDF with variables substituted with numbers

        """
        return self.Subs(expr).change_ring(CDF)

    def InitOperators(self):
        """
        Computes H, H_blocks, H_e, L by simple substitution

        """
        self.H = self.SubsNum(self.sym.H)
        self.H_blocks = [self.SubsNum(B) for B in self.sym.H_blocks]
        self.H_e = self.SubsNum(self.sym.H_e)
        self.L = self.sym.L.subs(gamma = self.params.gamma).change_ring(CDF)
        I2 = identity_matrix(CDF, 2)

    def RHS(self, rho):
        """
        Defines right-hand side of master equation ode

        """
        unitary_term = CDF(-I) * self.H.commutator(rho)

        L = self.L
        Lc = L.conjugate_transpose()
        lindblad_term = L * rho * Lc - 0.5 * rho.anticommutator(Lc * L)

        return unitary_term + lindblad_term

    def METimeStep(self, rho):
        """
        Computes time evolution of rho over time dt

        """
        rho = array(rho)
        matshape = rho.shape
        rho = rho.ravel()

        def f(t, y):
            """
            Defines right-hand side of master equation ode

            Can be used with scipy.integrate.ode ode solver

            """
            rho = matrix(y.reshape(matshape))
            rhs = self.RHS(rho)
            return array(rhs).ravel()

        dt = RDF(self.params.time_end / self.params.time_steps)

        r = ode(f)
        r.set_integrator('zvode')
        r.set_initial_value(rho)
        r.integrate(dt)

        rho = r.y.reshape(matshape)

        return matrix(rho)

    def ComputeTimeEvolution(self):
        """
        Computes time evolution

        """
        psi = full_basis_state(self.params.initial_state)
        rho = vec2dm(psi)

        self.rho_full_list = [rho]
        for t in range(1, self.params.time_steps):
            rho = self.METimeStep(rho)
            self.rho_full_list.append(rho)

        self.rho_list = [partial_trace_chain(x) for x in self.rho_full_list]
        self.rho_sink_list = [partial_trace_sink(x) for x in self.rho_full_list]

    def Rho(self, t):
        """
        Return density matrix at time t

        Evolution must have been computed beforehand

        """
        return self.rho_list[t]
