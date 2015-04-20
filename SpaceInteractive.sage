class SpaceInteractive(SpaceBase):
    """
    Makes interactive tables

    """
    def SubspaceSelector(self):
        result = selector(
            label = 'Subspace: ',
            values = range(self.qubits_count + 1),
            default = 1,
            buttons = True)
        return result

    def StateSlider(self, default = 1):
        return slider(self.states, default = default)

    def Bits(self, state):
        return bits(n = state, min_width = self.qubits_count)

    def BitsHeader(self):
        """
        Returns table header for bit string output

        """
        header = []
        for i in range(self.chain_len):
            header.append('ph{}'.format(i + 1))
            header.append('at{}'.format(i + 1))
        header.append('sink')

        return header

    def InteractiveState(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(state = self.StateSlider(),
                mode = ['transformed', 'full']):
            """
            Prints hilbert space basis vectors and excitation numbers
            User specifies initial state and basis

            Should be passed to interact()

            """
            choice = state
            if (mode == 'transformed'):
                state = self.T_rows[choice]
                state_e = choice
            else:
                state = choice
                state_e = self.T_columns[choice]

            subspace = self.GetSubspaceByState(state)

            html('<h2>Basis vector numbers</h2>')
            header = ['index', 'index_e', 'subspace', 'subspace_index',
                'subspace_size']
            row = [state, state_e, subspace.number, subspace.Index(state),
                subspace.Size()]
            html.table([row], header = header)

            html('<h2>Bits</h2>')
            row = self.Bits(state)
            html.table([row], header = self.BitsHeader())

        return inner

    def InteractiveSubspace(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(E = self.SubspaceSelector()):
            """
            Prints hilbert space subspace basis vectors
            User specifies subspace

            Should be passed to interact()

            """
            subspace = self.GetSubspaceByNumber(E)
            header = ['subspace_index', 'index_e', 'index']
            header.extend(self.BitsHeader())
            rows = []
            for i in subspace.states:
                row = [subspace.Index(i), self.T_columns[i], i]
                row.extend(self.Bits(i))
                rows.append(row)
            html.table(rows, header = header)

        return inner
