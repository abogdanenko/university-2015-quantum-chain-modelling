class SymbolicComputation:
    """
    Computes hamiltonian and other matrices symbolically
    """

    def ToEnergyBasis(self, A):
        """
        Returnes matrix A in energy basis
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

        self.H = self.H_sum.tensor_product(I4) \
            + I4.tensor_product(self.H_sum) \
            + self.H_tun
    
        self.H_blocks = []
        for E in energy_list:
            I = J = e_states(E)
            self.H_blocks.append(self.H[I, J])

        self.T = coordinate_change_matrix()

        self.H_e = self.ToEnergyBasis(self.H)

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
        l.append((r'H', self.H))
        
        l.append((r'T', self.T))
        l.append((r'H^e', self.H_e))        

        for E in energy_list:
            l.append((r'H_{}'.format(E), self.H_blocks[E]))

        html_vars(l)

    def ShowBasisStatesTableHTML(self):
        """
        Prints hilbert space basis vectors and their energies

        Works inside notebook interface
        """

        header = ['state', 'ph1', 'at1', 'ph2', 'at2', r'$E_{\rm total}$']
        header_e = header + [r'$\langle E \rangle$']
            
        html('<h2>Basis states</h2>')
        
        rows = []
        for i in states_list:
            row = [i]
            row.extend(bits(i, qubits_count))
            row.append(energy(i))
            row.append(self.H[i, i])
            rows.append(row)

        html.table(rows, header = header_e)
        
        rows = []
        for E in energy_list:
            html('<h2>Subspace (E = {})</h2>'.format(E))
            rows = []
            for i in e_states(E):
                row = [i]
                row.extend(bits(i, qubits_count))
                row.append(energy(i))
                rows.append(row)
            html.table(rows, header = header)

    def ComputeEigenVectors(self):
        """
        Computes eigenvalues and eigenvectors of each block
  
        Assumes omega_a = omega_c = omega
        """

        omega = SR.var('omega', domain = 'real')

        self.ev_blocks = []
        self.eigenvalues_blocks =[]
        for B in self.H_blocks:
            B1 = B.subs(omega_a = omega, omega_c = omega)
            ev = B1.eigenvectors_right()
            ev = flatten_eigen(ev)
            self.ev_blocks.append(ev)
            self.eigenvalues_blocks.append([value for value, vector in ev])

        self.ev_matrix = matrix(SR, states_count)
        column = 0
        row = 0
        for ev in self.ev_blocks:
            for value, vector in ev:
                for i in range(len(vector)):
                    self.ev_matrix[row + i, column] = vector[i]
                column += 1
            row += len(vector)


            
    def ShowEigenHTML(self):
        """
        Prints eigenvalues and eigenvectors of each block in HTML

        Works in notebook interface.
        """

        for E in energy_list:
            html(('<h3>Eigenvalues and eigenvectors of'
                ' $H_{}$: </h3>').format(E))
            show_eigen_html(self.ev_blocks[E])

