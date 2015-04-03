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
        Computes H, L by simple substitution

        """
        self.H = self.SubsNum(self.sym.H)
        self.L = self.SubsNum(self.sym.L)

    def ComputeTimeEvolution(self):
        """
        Computes time evolution

        """
        psi = basis_state(self.params.initial_state)
        rho = vec2dm(psi)
        E = exc_number(self.params.initial_state)

        integrator = MEIntegrator(
            rho = get_block(rho, E),
            H = get_block(self.H, E),
            L = get_block(self.L, E),
            dt = self.params.Dt())

        self.rho_subspace_list = integrator.Integrate(self.params.time_steps)

        self.rho_list = [get_full(x, E) for x in self.rho_subspace_list]
        self.rho_sink11 = [partial_trace_sink11(x, E) for x in self.rho_subspace_list]

    def Rho(self, t):
        """
        Return density matrix at time t

        Evolution must have been computed beforehand

        """
        return self.rho_list[t]
