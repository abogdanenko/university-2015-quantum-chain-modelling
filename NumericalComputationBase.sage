class NumericalComputationBase(ComputationBase):
    """
    Computes time evolution numerically

    """
    def __init__(self, space):
        super(NumericalComputationBase, self).__init__(space = space, ring = CDF)
        self.params = NumericalParams()
        self.ParamsChanged()

    def IterationTime(self, t):
        """
        Returns moment of time corresponding to iteration number t

        """
        return t / self.params.time_steps * self.params.time_end

    def ParamsChanged(self):
        """
        Computes operators and subspace

        """
        self.subspace = self.space.GetSubspaceByState(self.params.initial_state)

        self.ComputeOperators(
            omega_a = self.params.omega_a,
            omega_c = self.params.omega_c,
            alpha = self.params.alpha,
            beta = self.params.beta,
            gamma = self.params.gamma
        )

    def ComputeTimeEvolution(self):
        """
        Computes time evolution

        """
        rho_full = self.space.GetBasisDM(self.params.initial_state)

        integrator = MEIntegrator(
            rho = self.subspace.GetBlock(rho_full),
            H = self.subspace.GetBlock(self.H),
            L = self.subspace.GetBlock(self.L),
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

    def Conductivity(self):
        """
        Return chain conductivity

        """
        return RDF(average(array(self.rho_sink11)))
