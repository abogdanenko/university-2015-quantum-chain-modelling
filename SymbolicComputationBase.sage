class SymbolicComputationBase(ComputationBase):
    """
    Computes hamiltonian and lindblad operators symbolically

    """
    def __init__(self, space):
        super(SymbolicComputationBase, self).__init__(
            space = space,
            ring = SR,
            omega_a = SR.var('omega_a', domain = 'positive'),
            omega_c = SR.var('omega_c', domain = 'positive'),
            alpha = SR.var('alpha', domain = 'positive'),
            beta = SR.var('beta', domain = 'positive'),
            gamma = SR.var('gamma', domain = 'positive')
        )
