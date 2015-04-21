class SymbolicComputationInteractive(SymbolicComputationBase):
    """
    Makes interactive tables

    """
    def InteractiveSubspaceOperators(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(index = self.space.SubspaceSelector()):
            """
            Prints H, L for given subspace

            User specifies subspace

            Should be passed to interact()

            """
            subspace = self.space.subspaces[index]

            name = r'H_{}'.format(index)
            var = subspace.GetBlock(self.H)
            l = [(name, var)]

            name = r'L_{}'.format(index)
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

            (r'T', self.space.T),
            (r'H^{\rm ex}', self.H_e),

            (r'L', self.L)
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
