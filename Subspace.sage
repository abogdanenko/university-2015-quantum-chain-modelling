class Subspace(object):
    """
    Represents hilbert space subspace

    """
    selector = selector(
        label = 'Subspace: ',
        values = exc_list,
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
