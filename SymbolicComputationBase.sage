class SymbolicComputationBase(object):
    """
    Computes hamiltonian and lindblad operators symbolically

    """
    def __init__(self, space):
        self.space = space

        omega_a = SR.var('omega_a', domain = 'positive')
        omega_c = SR.var('omega_c', domain = 'positive')
        alpha = SR.var('alpha', domain = 'positive')
        beta = SR.var('beta', domain = 'positive')

        self.psi_g = matrix(SR, 2, 1, [1, 0])
        self.psi_e = matrix(SR, 2, 1, [0, 1])
        self.psi_ph0 = self.psi_g
        self.psi_ph1 = self.psi_e
        self.a = self.psi_ph0 * self.psi_ph1.conjugate_transpose()
        self.a_plus = self.a.conjugate_transpose()
        self.sigma_minus = self.a
        self.sigma_plus = self.sigma_minus.conjugate_transpose()
        self.H_field = omega_c * self.a_plus * self.a
        self.H_at = omega_a * self.sigma_plus * self.sigma_minus
        self.H_field_at = alpha * (self.a.tensor_product(self.sigma_plus) \
            + self.a_plus.tensor_product(self.sigma_minus))
        I2 = identity_matrix(SR, 2)

        self.H_tun = beta * (self.a.tensor_product(I2).tensor_product(
            self.a_plus).tensor_product(I2) + self.a_plus.tensor_product(
            I2).tensor_product(self.a).tensor_product(I2))

        self.H_sum = self.H_field.tensor_product(I2) \
            + I2.tensor_product(self.H_at) + self.H_field_at

        self.H_chain = matrix(SR, 2 ** (self.space.chain_len * 2))
        for i in range(self.space.chain_len):
            left = identity_matrix(SR, 2 ** (2 * i))
            right = identity_matrix(SR, 2 ** (2 * (self.space.chain_len - i - 1)))
            self.H_chain += left.tensor_product(self.H_sum).tensor_product(right)

        for i in range(self.space.chain_len - 1):
            left = identity_matrix(SR, 2 ** (2 * i))
            right = identity_matrix(SR, 2 ** (2 * (self.space.chain_len - i - 2)))
            self.H_chain += left.tensor_product(self.H_tun).tensor_product(right)

        self.H = self.H_chain.tensor_product(I2)
        self.H_e = self.space.ToExcBasis(self.H)

        gamma = SR.var('gamma', domain = 'positive')
        In1 = identity_matrix(SR, 2 ** (self.space.qubits_count - 2))

        self.L = gamma * In1.tensor_product(
            self.sigma_minus).tensor_product(self.sigma_plus)

        self.L_e = self.space.ToExcBasis(self.L)
