class SpaceCommon(object):
    """
    Common methods of SpaceBase and Subspace

    """
    def Size(self):
        """
        Returns number of states

        """
        return len(self.states)

    def StatesIndices(self):
        """
        Returns list of indices of states

        """
        return [state.index for state in self.states]

    def StateSlider(self, default = None):
        """
        Returns slider of states

        To be used inside interact

        """

        return slider(self.states, label = 'State:', default = default)
