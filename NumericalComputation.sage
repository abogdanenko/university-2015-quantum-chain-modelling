class NumericalComputation:
    """
    Computes time evolution numerically
    """

    def __init__(self, sym):
        self.sym = sym
        self.params = NumericalParams()

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

    def ComputeHamiltonian(self):
        """
        Computes H, H_blocks, H_e by simple substitution
        """ 

        self.H = self.SubsNum(self.sym.H)
        self.H_blocks = [self.SubsNum(B) for B in self.sym.H_blocks]
        self.H_e = self.SubsNum(self.sym.H_e)

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

    def ComputeTimeEvolution(self):
        """
        Computes time evolution matrix at each time t

        Also, takes squared module of each element of the time evolution matrix
        """

        self.U = [self.TimeEvolutionMatrix(self.IterationTime(t))
            for t in range(self.params.time_steps)]

        self.U_norm = []
        for m in self.U:
            n = m.nrows()
            r = matrix(RDF, n, n)
            for i in range(n):
                for j in range(n):
                    r[i, j] = norm(m[i, j])
            self.U_norm.append(r)

        self.U_e_norm = [self.sym.ToEnergyBasis(m) for m in self.U_norm]

    def DiagDist(self, initial_state, t):
        """
        Returns distance between state density matrix and a set of diagonal matrices

        State is a state at time step t. Distance between matrices is frobenius norm
        of matrix difference.
        """

        s = RDF()
        for i in states_list:
            for j in states_list:
                if (i != j):
                    U = self.U[t]
                    psi_i = U[i, initial_state]
                    psi_j = U[j, initial_state]
                    s += norm(psi_i * psi_j)
        return sqrt(s)

    def Entropy1(self, t, initial_state):
        """
        Returnes Von Neumann entropy of reduced density matrix
        """
        m = self.U[t]
        c = m.column(initial_state)
        y = PsiReduced1(c).column()
        rho = y * y.conjugate_transpose()
        # drop imag part, it should be zero
        ev = map(lambda z: z.real(), rho.eigenvalues())
        return sum([xlnx(p) for p in ev])
    
