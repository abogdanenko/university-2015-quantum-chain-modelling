class SymbolicComputationInteractive(SymbolicComputationBase):
    """
    Makes interactive tables

    """
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

            html('<h2>Basis vector numbers</h2>')
            header = ['index', 'index_e', 'subspace', 'subspace_index',
                'subspace_size']
            row = [state, state_e, exc_number(state), exc_index(state),
                exc_size(state)]
            html.table([row], header = header)

            html('<h2>Bits</h2>')
            header = ['ph1', 'at1', 'ph2', 'at2', 'sink']
            row = bits(state, qubits_count)
            html.table([row], header = header)

        return inner

    def InteractiveSubspace(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(subspace = slider(exc_list, default = 1)):
            """
            Prints hilbert space subspace basis vectors
            User specifies subspace

            Should be passed to interact()

            """
            header = ['subspace_index', 'index_e', 'index',
                'ph1', 'at1', 'ph2', 'at2', 'sink']
            rows = []
            for i in e_states(subspace):
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
        def inner(subspace = slider(exc_list, default = 1)):
            """
            Prints H, L for given subspace

            User specifies subspace

            Should be passed to interact()

            """
            name = r'H_{}'.format(subspace)
            var = get_block(self.H, subspace)
            l = [(name, var)]

            name = r'L_{}'.format(subspace)
            var = get_block(self.L, subspace)
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
