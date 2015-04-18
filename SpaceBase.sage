class SpaceBase(object):
    """
    Represents hilbert space

    """
    def __init__(self, chain_len):
        self.chain_len = chain_len
        self.qubits_count = 2 * self.chain_len + 1
        self.states_count = 2 ** self.qubits_count
        self.states = range(self.states_count)
        self.InitTransform()

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
        for E in range(self.qubits_count + 1):
            self.T_rows.extend(self.GetSubspaceStates(E))

        self.T_columns = [0] * self.states_count
        self.T = matrix(self.states_count)

        for j in self.states:
            i = self.T_rows[j]
            self.T_columns[i] = j
            self.T[i, j] = 1

    def ToExcBasis(self, A):
        """
        Returnes matrix A in excitation basis

        """
        return self.T.transpose() * A * self.T

    def ExcCount(self, state):
        """
        Returns total number of excitations

        """
        return self.Bits(state).count(1)

    def GetSubspaceStates(self, E):
        """
        Returns a list of states with total number of excitations E

        """
        return [i for i in self.states if self.ExcCount(i) == E]

    def GetSubspaceByNumber(self, number):
        """
        Returns subspace by number

        """
        return Subspace(number, self)

    def GetSubspaceByState(self, state):
        """
        Returns subspace to which state belongs

        """
        return self.GetSubspaceByNumber(self.ExcCount(state))

    def GetStateSubspaceIndex(state):
        """
        Returns index of state within its subspace

        """
        subspace = self.GetSubspaceByState(state)
        return subspace.Index(state)

    def GetBasisVector(self, state):
        """
        Returns basis state number ``state``

        """
        psi = vector(CDF, self.states_count)
        psi[state] = 1
        return psi

    def GetBasisDM(self, state):
        """
        Returns basis state number ``state`` as density matrix

        """
        v = self.GetBasisVector(state).column()
        return v * v.conjugate_transpose()
