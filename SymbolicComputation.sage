class SymbolicComputation(object):
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

        I4 = identity_matrix(SR, 4)

        self.H_chain = self.H_sum.tensor_product(I4) \
            + I4.tensor_product(self.H_sum) \
            + self.H_tun

        self.H = self.H_chain.tensor_product(I2)

        self.H_blocks = []
        for E in exc_list:
            I = J = e_states(E)
            self.H_blocks.append(self.H[I, J])

        self.T = coordinate_change_matrix()

        self.H_e = self.ToExcBasis(self.H)

        gamma = SR.var('gamma', domain = 'positive')
        self.L = gamma * I4.tensor_product(I2).tensor_product(
            self.sigma_minus).tensor_product(self.sigma_plus)

    def ShowVarsHTML(self):
        """
        Prints matrices of quantum-mech. operators used in computations

        Works inside notebook interface

        """
        l = []

        l.append((r'|g\rangle', self.psi_g))
        l.append((r'|e\rangle', self.psi_e))
        l.append((r'|0\rangle', self.psi_ph0))
        l.append((r'|1\rangle', self.psi_ph1))

        l.append((r'a', self.a))
        l.append((r'a^{+}', self.a_plus))
        l.append((r'\sigma^{-}', self.sigma_minus))
        l.append((r'\sigma^{+}', self.sigma_plus))

        l.append((r'H_{\rm field}', self.H_field))
        l.append((r'H_{\rm at}', self.H_at))
        l.append((r'H_{\rm field,at}', self.H_field_at))
        l.append((r'H_{\rm sum}', self.H_sum))
        l.append((r'H_{\rm tun}', self.H_tun))
        l.append((r'H_{\rm chain}', self.H_chain))
        l.append((r'H', self.H))

        l.append((r'T', self.T))
        l.append((r'H^{\rm ex}', self.H_e))

        for E in exc_list:
            l.append((r'H_{}'.format(E), self.H_blocks[E]))

        l.append((r'L', self.L))

        html_vars(l)

    def ShowBasisStatesTableHTML(self):
        """
        Prints hilbert space basis vectors, excitation numbers and energies

        Works inside notebook interface

        """
        header = ['state', 'ph1', 'at1', 'ph2', 'at2', r'$N_{\rm ex}$']
        header_e = header + [r'$\langle H \rangle$']

        html('<h2>Basis states</h2>')

        rows = []
        for i in states_list:
            row = [i]
            row.extend(bits(i, qubits_count))
            row.append(exc_number(i))
            row.append(self.H[i, i])
            rows.append(row)

        html.table(rows, header = header_e)

        rows = []
        for E in exc_list:
            html(r'<h2>Subspace ($N_{{\rm ex}} = {}$)</h2>'.format(E))
            rows = []
            for i in e_states(E):
                row = [i]
                row.extend(bits(i, qubits_count))
                row.append(exc_number(i))
                rows.append(row)
            html.table(rows, header = header)