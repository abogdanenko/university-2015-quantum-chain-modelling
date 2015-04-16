class Subspace(object):
    """
    Represents hilbert space subspace

    """
    selector = selector(
        label = 'Subspace: ',
        values = range(qubits_count + 1),
        default = 1,
        buttons = True)

    def __init__(self, number):
        self.number = number
        self.states = e_states(number)
        self.size = len(self.states)
        self.colors = rainbow(self.size)

    def GetBlock(self, A):
        """
        Returns subspace block of A in exc basis

        """
        return get_block(A, self.number)

    def GetFull(self, B):
        """
        Returns full matrix A given subspace block of A

        Sets all elements of A to zero except elements in block B

        """
        A = matrix(B.base_ring(), states_count)
        I = J = self.states
        A[I, J] = B
        return A

    def partial_trace_sink11(self, rho):
        """
        Returns rho_sink[1,1]

        """
        s = 0
        for state in self.states:
            if state % 2: # sink bit is set
                i = exc_index(state)
                s += rho[i, i]
        return abs(s)
