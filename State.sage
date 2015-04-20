class State(object):
    """
    Represents hilbert space basis state

    """
    def __init__(self, index):
        self.index = index

    def SetSpace(self, space):
        """
        Sets space

        """
        self.space = space

    def SetSubspace(self, subspace):
        """
        Sets subspace

        """
        self.subspace = subspace

    def ExcCount(self):
        """
        Returns total number of excitations

        """
        return bits(self.index).count(1)
