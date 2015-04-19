class NumericalComputationAnimation(NumericalComputationPlots):
    """
    Saves animated plots as files

    """
    def SaveRhoAnimation(self, mode = 'subspace'):
        """
        Saves density matrix plot png image for each moment of time

        """
        dirname = tmp_dir()

        for t in range(self.params.time_steps):
            filename = os.path.join(dirname, '{:08d}.png'.format(t))
            plot = self.PlotRho(t, mode)
            plot.save_image(filename)

        return dirname

    def SaveProbabilityBarChartAnimation(self, mode = 'subspace'):
        """
        Saves bar chart png image for each moment of time

        """
        dirname = tmp_dir()

        for t in range(self.params.time_steps):
            filename = os.path.join(dirname, '{:08d}.png'.format(t))
            self.ProbabilityBarChart(t, mode, filename)

        return dirname
