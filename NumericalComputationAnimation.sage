class NumericalComputationAnimation(object):
    """
    Saves animated plots as files for class NumericalComputation
    """

    def SaveUnitaryEvolutionMatrixAnimation(self, filename, basis_e = False):
        """
        Saves matrix plots of U as gif animation
        """

        l = [self.PlotUnitaryEvolutionMatrix(t, basis_e)
            for t in range(self.params.time_steps)]
        animation = animate(l)
        animation.gif(savefile = filename, show_path = True)
    
    def SaveRhoAnimation(self, filename, initial_state, basis_e = False):
        """
        Saves density matrix plot as gif animation
        """

        l = [self.PlotRho(initial_state, t, basis_e)
            for t in range(self.params.time_steps)]
        animation = animate(l)
        animation.gif(savefile = filename, show_path = True)
    
    def SaveProbabilityBarChartAnimation(self, dirname, initial_state, basis_e = False):
        """
        Saves bar chart png image for each moment of time
        """

        for t in range(self.params.time_steps):
            filename = os.path.join(dirname, '{:06d}.png'.format(t))
            self.ProbabilityBarChart(initial_state, t, basis_e, filename)
