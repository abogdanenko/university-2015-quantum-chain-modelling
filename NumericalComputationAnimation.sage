class NumericalComputationAnimation(NumericalComputationPlots):
    """
    Saves animated plots as files

    """
    def SaveUnitaryEvolutionMatrixAnimation(self, basis_e = False):
        """
        Saves U matrix plot png image for each moment of time

        """
        l = [self.PlotUnitaryEvolutionMatrix(t, basis_e)
            for t in range(self.params.time_steps)]
        animation = animate(l)
        return animation.png()

    def SaveRhoAnimation(self, basis_e = False):
        """
        Saves density matrix plot png image for each moment of time

        """
        l = [self.PlotRho(t, basis_e)
            for t in range(self.params.time_steps)]
        animation = animate(l)
        return animation.png()

    def SaveProbabilityBarChartAnimation(self, basis_e = False):
        """
        Saves bar chart png image for each moment of time

        """
        dirname = tmp_dir()

        for t in range(self.params.time_steps):
            filename = os.path.join(dirname, '{:08d}.png'.format(t))
            self.ProbabilityBarChart(t, basis_e, filename)

        return dirname
