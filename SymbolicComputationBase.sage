class SymbolicComputationBase(object):
    """
    Computes hamiltonian symbolically

    """
    def ToExcBasis(self, A):
        """
        Returnes matrix A in excitation basis

        """
        return self.T.transpose() * A * self.T

    def __init__(self):
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

        self.H_chain = matrix(SR, 2 ** (chain_len * 2))
        for i in range(chain_len):
            left = identity_matrix(SR, 2 ** (2 * i))
            right = identity_matrix(SR, 2 ** (2 * (chain_len - i - 1)))
            self.H_chain += left.tensor_product(self.H_sum).tensor_product(right)

        for i in range(chain_len - 1):
            left = identity_matrix(SR, 2 ** (2 * i))
            right = identity_matrix(SR, 2 ** (2 * (chain_len - i - 2)))
            self.H_chain += left.tensor_product(self.H_tun).tensor_product(right)

        self.H = self.H_chain.tensor_product(I2)
        self.InitTransform()
        self.H_e = self.ToExcBasis(self.H)

        gamma = SR.var('gamma', domain = 'positive')
        In1 = identity_matrix(SR, 2 ** (qubits_count - 2))

        self.L = gamma * In1.tensor_product(
            self.sigma_minus).tensor_product(self.sigma_plus)

        self.L_e = self.ToExcBasis(self.L)

    def InitTransform(self):
        """
        Initializes T, T_rows, T_columns

        T - coordinate change matrix to excitation basis

        Invariants:
        T[i, T_columns[i]] = 1
        T[T_rows[j], j] = 1
        T_rows[T_columns[i]] = i
        T_columns[T_rows[j]] = j

        """
        self.T_rows = []
        for E in exc_list:
            self.T_rows.extend(e_states(E))

        self.T_columns = [0] * states_count
        self.T = matrix(states_count)

        for j in states_list:
            i = self.T_rows[j]
            self.T_columns[i] = j
            self.T[i, j] = 1
