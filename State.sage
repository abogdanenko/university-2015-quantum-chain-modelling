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

    def SubspaceIndex(self):
        """
        Returns index of state within its subspace

        """
        return self.subspace.Index(self)

    def Vector(self):
        """
        Returns vector representation of the state

        """
        psi = vector(CDF, self.space.StatesCount())
        psi[self.Index()] = 1
        return psi

    def DensityMatrix(self):
        """
        Returns density matrix representation of the state

        """
        v = self.Vector().column()
        return v * v.conjugate_transpose()

    def Bits(self):
        """
        Returns list of bits

        """
        return bits(n = self.index, min_width = self.space.QubitsCount())
