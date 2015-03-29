class NumericalComputationAnimation(NumericalComputationPlots):
    """
    Saves animated plots as files

    """
    def SaveRhoAnimation(self, basis_e = False):
        """
        Saves density matrix plot png image for each moment of time

        """
        dirname = tmp_dir()

        for t in range(self.params.time_steps):
            filename = os.path.join(dirname, '{:08d}.png'.format(t))
            plot = self.PlotRho(t, basis_e)
            plot.save_image(filename)

        return dirname

    def SaveRhoFullAnimation(self, basis_e = False):
        """
        Saves density matrix plot png image for each moment of time

        """
        dirname = tmp_dir()

        for t in range(self.params.time_steps):
            filename = os.path.join(dirname, '{:08d}.png'.format(t))
            plot = self.PlotRhoFull(t, basis_e)
            plot.save_image(filename)

        return dirname

    def SaveProbabilityBarChartAnimation(self, basis_e = False):
        """
        Saves bar chart png image for each moment of time

        """
        dirname = tmp_dir()

        for t in range(self.params.time_steps):
            filename = os.path.join(dirname, '{:08d}.png'.format(t))
            self.ProbabilityBarChart(t, basis_e, filename)

        return dirname
