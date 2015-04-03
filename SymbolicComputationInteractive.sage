class SymbolicComputationInteractive(SymbolicComputationBase):
    """
    Makes interactive tables

    """
    def BitsHeader(self):
        """
        Returns table header for bit string output

        """
        header = []
        for i in range(chain_len):
            header.append('ph{}'.format(i + 1))
            header.append('at{}'.format(i + 1))
        header.append('sink')

        return header

    def InteractiveState(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(state = slider(states_list),
                basis_e = [True, False]):
            """
            Prints hilbert space basis vectors and excitation numbers
            User specifies initial state and basis

            Should be passed to interact()

            """
            choice = state
            if basis_e:
                state = self.T_rows[choice]
                state_e = choice
            else:
                state = choice
                state_e = self.T_columns[choice]

            subspace = get_subspace(state)

            html('<h2>Basis vector numbers</h2>')
            header = ['index', 'index_e', 'subspace', 'subspace_index',
                'subspace_size']
            row = [state, state_e, subspace.number, exc_index(state),
                subspace.size]
            html.table([row], header = header)

            html('<h2>Bits</h2>')
            row = bits(state, qubits_count)
            html.table([row], header = self.BitsHeader())

        return inner

    def InteractiveSubspace(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(E = Subspace.selector):
            """
            Prints hilbert space subspace basis vectors
            User specifies subspace

            Should be passed to interact()

            """
            subspace = Subspace(E)
            header = ['subspace_index', 'index_e', 'index']
            header.extend(self.BitsHeader())
            rows = []
            for i in subspace.states:
                row = [exc_index(i), self.T_columns[i], i]
                row.extend(bits(i, qubits_count))
                rows.append(row)
            html.table(rows, header = header)

        return inner

    def InteractiveSubspaceOperators(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(E = Subspace.selector):
            """
            Prints H, L for given subspace

            User specifies subspace

            Should be passed to interact()

            """
            subspace = Subspace(E)

            name = r'H_{}'.format(subspace.number)
            var = subspace.GetBlock(self.H)
            l = [(name, var)]

            name = r'L_{}'.format(subspace.number)
            var = subspace.GetBlock(self.L)
            l.append((name, var))

            html_vars(l)

        return inner

    def InteractiveVars(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        l = [
            (r'|g\rangle', self.psi_g),
            (r'|e\rangle', self.psi_e),
            (r'|0\rangle', self.psi_ph0),
            (r'|1\rangle', self.psi_ph1),

            (r'a', self.a),
            (r'a^{+}', self.a_plus),
            (r'\sigma^{-}', self.sigma_minus),
            (r'\sigma^{+}', self.sigma_plus),

            (r'H_{\rm field}', self.H_field),
            (r'H_{\rm at}', self.H_at),
            (r'H_{\rm field,at}', self.H_field_at),
            (r'H_{\rm sum}', self.H_sum),
            (r'H_{\rm tun}', self.H_tun),
            (r'H_{\rm chain}', self.H_chain),
            (r'H', self.H),

            (r'T', self.T),
            (r'H^{\rm ex}', self.H_e),

            (r'L', self.L),

            (r'L^{\rm ex}', self.L_e)
        ]

        d = dict(l)

        latex_vars = ['${}$'.format(x[0]) for x in l]

        variable_selector = selector(
            label = 'Show variable: ',
            values = latex_vars,
            buttons = True,
            ncols = 4,
            nrows = 5)

        def inner(variable = variable_selector):
            """
            Prints matrix chosen by the user

            Should be passed to interact()

            """
            key = variable[1:-1]
            value = d[key]
            html_var(key, value)

        return inner
