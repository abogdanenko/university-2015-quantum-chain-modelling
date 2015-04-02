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
            omega_c = self.params.omega_c,
            gamma = self.params.gamma)
        return result

    def SubsNum(self, expr):
        """
        Returns expr coerced to CDF with variables substituted with numbers

        """
        return self.Subs(expr).change_ring(CDF)

    def InitOperators(self):
        """
        Computes H, H_blocks, L, L_blocks by simple substitution

        """
        self.H = self.SubsNum(self.sym.H)
        self.H_blocks = [self.SubsNum(B) for B in self.sym.H_blocks]
        self.L = self.SubsNum(self.sym.L)
        self.L_blocks = [self.SubsNum(B) for B in self.sym.L_blocks]

    def ComputeTimeEvolution(self):
        """
        Computes time evolution

        """
        psi = basis_state(self.params.initial_state)
        rho = vec2dm(psi)

        dt = RDF(self.params.time_end / self.params.time_steps)
        integrator = MEIntegrator(rho = rho, H = self.H, L = self.L, dt = dt)

        self.rho_list = integrator.Integrate(self.params.time_steps)

        self.rho_sink_list = [partial_trace_sink(x) for x in self.rho_list]

    def Rho(self, t):
        """
        Return density matrix at time t

        Evolution must have been computed beforehand

        """
        return self.rho_list[t]
