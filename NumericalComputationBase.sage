from scipy.integrate import ode
from numpy import array

class NumericalComputationBase(object):
    """
    Computes time evolution numerically

    """
    def __init__(self, sym):
        self.sym = sym
        self.params = NumericalParams(self.sym.unitary)
        self.InitOperators()
        self.ComputeEigenVectors()

    def TimeEvolutionMatrix(self, t):
        """
        Returns time evolution matrix at time t

        """
        return exp(CDF(-I * t) * self.H)

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

        if not self.sym.unitary:
            self.L = self.sym.L.subs(gamma = self.params.gamma).change_ring(CDF)

    def ComputeEigenVectors(self):
        """
        Computes eigenvalues and eigenvectors of each block

        """
        self.ev_blocks = [ev_flat_sorted(B) for B in self.H_blocks]

    def ShowEigenHTML(self):
        """
        Prints eigenvalues and eigenvectors of each block in HTML

        Works in notebook interface.

        """
        for E in energy_list:
            html(('<h3>Eigenvalues and eigenvectors of $H_{}$: </h3>'
                ).format(E))
            show_eigen_html(self.ev_blocks[E])


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

    def ComputeTimeEvolution(self, initial_states = []):
        """
        Computes time evolution

        Accepts a set of initial states to compute evolution for in master
        equation case

        """
        if self.sym.unitary:
            self.U = [self.TimeEvolutionMatrix(self.IterationTime(t))
                for t in range(self.params.time_steps)]

        else:
            rho_initial = [vec2dm(basis_state(state)) for state in states_list]

            self.rho_list = [rho_initial]
            for t in range(1, self.params.time_steps):
                l = []
                for state in states_list:
                    rho = self.rho_list[-1][state]
                    if state in initial_states:
                        rho = self.METimeStep(rho)
                    l.append(rho)

                self.rho_list.append(l)


    def Rho(self, initial_state, t):
        """
        Return density matrix at time t

        Evolution must have been computed beforehand

        """
        if self.sym.unitary:
            U = self.U[t]
            vec = U.column(initial_state)
            psi = vec.column()
            return psi * psi.conjugate_transpose()
        else:
            return self.rho_list[t][initial_state]

    def DiagDist(self, initial_state, t):
        """
        Returns distance between density matrix and a set of diagonal matrices

        State is a state at time step t. Distance between matrices is frobenius
        norm of matrix difference.

        """
        s = RDF()
        rho = self.Rho(initial_state, t)
        for i in states_list:
            for j in states_list:
                if (i != j):
                    s += norm(rho[i, j])
        return sqrt(s)

    def Entropy1(self, initial_state, t):
        """
        Returnes Von Neumann entropy of reduced density matrix

        """
        rho = self.Rho(initial_state, t)
        rho1 = partial_trace1(rho)
        # drop imag part, it should be zero
        ev = map(abs, rho1.eigenvalues())
        return -1 * sum([xlnx(p) for p in ev])
