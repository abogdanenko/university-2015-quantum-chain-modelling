class SymbolicComputationBase(ComputationBase):
    """
    Computes hamiltonian and lindblad operators symbolically

    """
    def __init__(self, space):
        super(SymbolicComputationBase, self).__init__(space = space, ring = SR)

        self.ComputeOperators(
            omega_a = SR.var('omega_a', domain = 'positive'),
            omega_c = SR.var('omega_c', domain = 'positive'),
            alpha = SR.var('alpha', domain = 'positive'),
            beta = SR.var('beta', domain = 'positive'),
            gamma_s = SR.var('gamma_s', domain = 'positive'),
            gamma_d = SR.var('gamma_d', domain = 'positive')
        )
