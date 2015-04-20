class SpaceBase(object):
    """
    Represents hilbert space

    """
    def __init__(self, chain_len):
        self.chain_len = chain_len
        self.InitStatesSubspaces()
        self.InitTransform()

    def InitStatesSubspaces(self):
        """
        Initializes states, subspaces

        """
        states_count = 2 ** self.QubitsCount()
        self.states = []
        subspaces = {}
        for i in range(states_count):
            state = State(i)
            subspace_index = state.ExcCount()
            if subspace_index not in subspaces:
                subspaces[subspace_index] = Subspace(subspace_index, self)
            subspaces[subspace_index].AddState(state)

            state.SetSpace(self)
            self.states.append(state)

        # get list of values ordered by key
        self.subspaces = [subspaces[i] for i in range(len(subspaces))]

    def QubitsCount(self):
        """
        Returns number of qubits

        """
        return 2 * self.chain_len + 1

    def StatesCount(self):
        """
        Returns number of states

        """
        return len(self.states)

    def SubspacesCount(self):
        """
        Returns number of subspaces

        """
        return len(self.subspaces)

    def StatesIndices(self):
        """
        Returns list of indices of states

        """
        return [state.index for state in self.states]

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
        for subspace in self.subspaces:
            self.T_rows.extend(subspace.StatesIndices())

        self.T_columns = [0] * self.StatesCount()
        self.T = matrix(self.StatesCount())

        for j in self.StatesIndices():
            i = self.T_rows[j]
            self.T_columns[i] = j
            self.T[i, j] = 1

    def ToExcBasis(self, A):
        """
        Returnes matrix A in excitation basis

        """
        return self.T.transpose() * A * self.T

    def GetSubspaceByIndex(self, index):
        """
        Returns subspace by index

        """
        return Subspace(index, self)

    def GetSubspaceByState(self, state):
        """
        Returns subspace to which state belongs

        """
        return self.GetSubspaceByIndex(self.ExcCount(state))
