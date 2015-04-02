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
