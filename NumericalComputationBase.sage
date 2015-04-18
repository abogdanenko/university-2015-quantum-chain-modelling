class NumericalComputationBase(object):
    """
    Computes time evolution numerically

    """
    def __init__(self, sym):
        self.sym = sym
        self.space = self.sym.space
        self.params = NumericalParams()
        self.ParamsChanged()

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

    def ParamsChanged(self):
        """
        Computes H, L by simple substitution

        """
        self.subspace = self.space.GetSubspace(self.params.initial_state)
        H = self.subspace.GetBlock(self.sym.H)
        L = self.subspace.GetBlock(self.sym.L)
        self.H = self.SubsNum(H)
        self.L = self.SubsNum(L)

    def ComputeTimeEvolution(self):
        """
        Computes time evolution

        """
        rho_full = self.space.GetBasisDM(self.params.initial_state)
        rho = self.subspace.GetBlock(rho_full)
        integrator = MEIntegrator(
            rho = rho,
            H = self.H,
            L = self.L,
            dt = self.params.Dt())

        self.rho = integrator.Integrate(self.params.time_steps)

        self.rho_sink11 = [self.subspace.partial_trace_sink11(rho)
            for rho in self.rho]

    def Rho(self, t):
        """
        Return density matrix at time t

        Evolution must have been computed beforehand

        """
        return self.subspace.GetFull(self.rho[t])
