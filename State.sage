class State(object):
    """
    Represents hilbert space basis state

    """
    def __init__(self, index):
        self.index = index

    def __repr__(self):
        return "ind={}, ind_e={}, ss={}, ind_ss={}".format(self.index,
            self.TransformedIndex(), self.subspace.index, self.SubspaceIndex())

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

    def TransformedIndex(self):
        """
        Returns index of state in e_basis

        """
        return self.space.T_columns[self.index]

    def Vector(self):
        """
        Returns vector representation of the state

        """
        psi = vector(CDF, self.space.Size())
        psi[self.index] = 1
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

    def Show(self):
        """
        Shows indices in default and transformed bases, shows bits

        """
        html('<h2>State indices</h2>')
        header = ['index', 'index_e', 'index_subspace', 'subspace',
            'subspace_size']
        row = [self.index, self.TransformedIndex(), self.SubspaceIndex(),
            self.subspace.index, self.subspace.Size()]
        html.table([row], header = header)

        html('<h2>Bits</h2>')
        html.table([self.Bits()], header = self.space.BitsHeader())
