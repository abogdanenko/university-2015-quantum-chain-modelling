class SpaceInteractive(SpaceBase):
    """
    Makes interactive tables

    """
    def SubspaceSelector(self):
        """
        Returns selector of subspaces

        To be used inside interact

        """
        result = selector(
            label = 'Subspace: ',
            values = range(self.SubspacesCount()),
            default = 1,
            buttons = True)
        return result

    def StateSlider(self, default = 1):
        return slider(self.states, default = default)

    def BitsHeader(self):
        """
        Returns table header for bit string output

        """
        header = []
        for i in range(self.chain_len):
            header.append('ph{}'.format(i + 1))
            header.append('at{}'.format(i + 1))
        header.append('sink')

        return header

    def InteractiveState(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(index = self.StateSlider()):
            """
            Shows state interactively

            Should be passed to interact()

            """
            self.states[index].Show()

        return inner

    def InteractiveSubspace(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(index = self.SubspaceSelector()):
            """
            Prints hilbert space subspace basis vectors
            User specifies subspace

            Should be passed to interact()

            """
            self.subspaces[index].Show()

        return inner
