class Subspace(SpaceCommon):
    """
    Represents hilbert space subspace

    """
    def __init__(self, index, space):
        self.index = index
        self.space = space
        self.states = []

    def AddState(self, state):
        """
        Adds state to list of states

        """
        state.SetSubspace(self)
        self.states.append(state)

    def Colors(self):
        """
        Returns list of colors, one color for one state

        """
        return rainbow(self.Size())

    def GetBlock(self, A):
        """
        Returns subspace block of A in exc basis

        """
        I = J = self.StatesIndices()
        return A[I, J]

    def GetFull(self, B):
        """
        Returns full matrix A given subspace block of A

        Sets all elements of A to zero except elements in block B

        """
        A = matrix(B.base_ring(), self.space.Size())
        I = J = self.StatesIndices()
        A[I, J] = B
        return A

    def Index(self, state):
        """
        Returns index of state within its subspace

        """
        return self.states.index(state)

    def partial_trace_sink11(self, rho):
        """
        Returns rho_sink[1,1]

        """
        s = 0
        for state in self.states:
            if state.index % 2: # sink bit is set
                i = state.SubspaceIndex()
                s += rho[i, i]
        return abs(s)

    def Show(self):
        """
        Shows states

        """
        header = ['index_subspace', 'index_e', 'index']
        header.extend(self.space.BitsHeader())
        rows = []
        for state in self.states:
            row = [state.SubspaceIndex(), state.TransformedIndex(), state.index]
            row.extend(state.Bits())
            rows.append(row)
        html.table(rows, header = header)
