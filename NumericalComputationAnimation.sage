class NumericalComputationAnimation(object):
    """
    Saves animated plots as files for class NumericalComputation

    """
    def SaveUnitaryEvolutionMatrixAnimation(self, basis_e = False):
        """
        Saves U matrix plot png image for each moment of time

        """
        l = [self.PlotUnitaryEvolutionMatrix(t, basis_e)
            for t in range(self.params.time_steps)]
        animation = animate(l)
        return animation.png()

    def SaveRhoAnimation(self, initial_state, basis_e = False):
        """
        Saves density matrix plot png image for each moment of time

        """
        l = [self.PlotRho(initial_state, t, basis_e)
            for t in range(self.params.time_steps)]
        animation = animate(l)
        return animation.png()

    def SaveProbabilityBarChartAnimation(self, initial_state, basis_e = False):
        """
        Saves bar chart png image for each moment of time

        """
        dirname = tmp_dir()

        for t in range(self.params.time_steps):
            filename = os.path.join(dirname, '{:08d}.png'.format(t))
            self.ProbabilityBarChart(initial_state, t, basis_e, filename)

        return dirname
