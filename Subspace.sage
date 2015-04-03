class Subspace(object):
    """
    Represents hilbert space subspace

    """
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
        Returns full matrix A given E-th block of A in exc basis

        Sets all elements of A to zero except elements in block B

        """
        return get_full(B, self.number)
