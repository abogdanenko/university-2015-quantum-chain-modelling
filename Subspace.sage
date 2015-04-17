class Subspace(object):
    """
    Represents hilbert space subspace

    """
    def __init__(self, number, space):
        self.number = number
        self.space = space
        self.states = self.space.GetSubspaceStates(number)
        self.size = len(self.states)
        self.colors = rainbow(self.size)

    def GetBlock(self, A):
        """
        Returns subspace block of A in exc basis

        """
        I = J = self.states
        return A[I, J]

    def GetFull(self, B):
        """
        Returns full matrix A given subspace block of A

        Sets all elements of A to zero except elements in block B

        """
        A = matrix(B.base_ring(), states_count)
        I = J = self.states
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
            if state % 2: # sink bit is set
                i = self.Index(state)
                s += rho[i, i]
        return abs(s)
