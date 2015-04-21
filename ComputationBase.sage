class ComputationBase(object):
    """
    Computes hamiltonian and lindblad operators given ring and parameters

    """
    def __init__(self, space, ring):
        self.space = space
        self.ring = ring
        self.InitBasisVectors()
        self.InitOperatorsASigma()

    def InitBasisVectors(self):
        """
        Initializes atom and photon states for one cavity

        """
        self.psi_g = matrix(self.ring, 2, 1, [1, 0])
        self.psi_e = matrix(self.ring, 2, 1, [0, 1])
        self.psi_ph0 = self.psi_g
        self.psi_ph1 = self.psi_e

    def InitOperatorsASigma(self):
        """
        Initializes operators a, a_plus, sigma_minus, sigma_plus

        """
        self.a = self.psi_ph0 * self.psi_ph1.conjugate_transpose()
        self.a_plus = self.a.conjugate_transpose()
        self.sigma_minus = self.a
        self.sigma_plus = self.sigma_minus.conjugate_transpose()

    def SandwichOperator(self, operator, qubits_left, qubits_right):
        """
        Returns operator acting in larger space

        """
        left = identity_matrix(self.ring, 2 ** qubits_left)
        right = identity_matrix(self.ring, 2 ** qubits_right)
        return left.tensor_product(operator).tensor_product(right)

    def ComputeHamiltonian(self, omega_a, omega_c, alpha, beta):
        """
        Computes hamiltonian of one cavity and of the whole system

        """
        self.H_field = omega_c * self.a_plus * self.a
        self.H_at = omega_a * self.sigma_plus * self.sigma_minus
        self.H_field_at = alpha * (self.a.tensor_product(self.sigma_plus) \
            + self.a_plus.tensor_product(self.sigma_minus))
        I2 = identity_matrix(self.ring, 2)

        self.H_tun = beta * (self.a.tensor_product(I2).tensor_product(
            self.a_plus).tensor_product(I2) + self.a_plus.tensor_product(
            I2).tensor_product(self.a).tensor_product(I2))

        self.H_sum = self.H_field.tensor_product(I2) \
            + I2.tensor_product(self.H_at) + self.H_field_at

        chain_qubits = self.space.QubitsCount() - 1
        self.H_chain = matrix(self.ring, 2 ** chain_qubits)

        for i in range(self.space.chain_len):
            self.H_chain += self.SandwichOperator(self.H_sum, 2 * i,
                chain_qubits - 2 * i - 2)

        for i in range(self.space.chain_len - 1):
            self.H_chain += self.SandwichOperator(self.H_tun, 2 * i,
                chain_qubits - 2 * i - 4)

        self.H = self.H_chain.tensor_product(I2)
        self.H_e = self.space.ToExcBasis(self.H)

    def ComputeLindbladOperators(self, gamma_s, gamma_d):
        """
        Computes lindblad operators

        """
        In1 = identity_matrix(self.ring, 2 ** (self.space.QubitsCount() - 2))

        L_sink = gamma_s * In1.tensor_product(
            self.sigma_minus).tensor_product(self.sigma_plus)

        self.L = [L_sink]

        for i in range(self.space.chain_len):
            L = self.SandwichOperator(
                operator = self.sigma_plus * self.sigma_minus,
                qubits_left = 2 * i + 1,
                qubits_right = self.space.QubitsCount() - 2 * i - 2)

            self.L.append(gamma_d * L)

    def ComputeOperators(self, omega_a, omega_c, alpha, beta, gamma_s, gamma_d):
        """
        Computes hamiltonian and lindblad operators

        """
        self.ComputeHamiltonian(omega_a, omega_c, alpha, beta)
        self.ComputeLindbladOperators(gamma_s, gamma_d)
