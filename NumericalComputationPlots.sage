import matplotlib.pyplot as plt

class NumericalComputationPlots(NumericalComputationBase):
    """
    Makes plots

    """
    def TimeList(self):
        """
        Returns list of time points

        """
        return [self.IterationTime(t) for t in range(self.params.time_steps)]

    def PlotState(self, index, color = 'red'):
        """
        Returns line plot of the state

        """
        state = self.space.states[index]
        i = state.SubspaceIndex()

        Y = []
        for t in range(self.params.time_steps):
            rho = self.rho[t]
            y = abs(rho[i, i])
            Y.append(y)

        l = zip(self.TimeList(), Y)

        plot_object = line(l,
            ymin = 0,
            ymax = 1,
            color = color,
            tick_formatter = 'latex',
            axes_labels = ['$t$', '$P$'],
            legend_label = r'$\rho_{{{0},{0}}}(t)$'.format(index))

        plot_object.set_legend_options(back_color = 'white', ncol = 2)

        return plot_object

    def PlotStates(self):
        """
        Returns multi-line plot of state evolution

        """
        pairs = zip(self.subspace.StatesIndices(), self.subspace.Colors())
        plots = [self.PlotState(*x) for x in pairs]
        return sum(plots)

    def ShowStates(self):
        """
        Shows multi-line plot of state evolution

        """
        html('<h2>Density matrix diagonal elems</h2>')
        show(self.PlotStates())

    def ProbabilityBarChart(
            self,
            t,
            mode = 'subspace',
            filename = 'ProbabilityBarChart'):
        """
        Returns bar chart of state at time t

        """
        rho = self.rho[t]
        x = range(self.subspace.Size())

        if (mode == 'transformed'):
            rho = self.space.ToExcBasis(self.Rho(t))
            x = self.space.StatesIndices()

        if (mode == 'full'):
            rho = self.Rho(t)
            x = self.space.StatesIndices()

        d = map(abs, rho.diagonal())
        ylabel = r'$\rho_{j,j}$'

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        title1 = 't = {:7.2f}'.format(float(num.IterationTime(t)))
        title2 = r'\rho(0) = |{0}\rangle \langle {0}|'.format(
                 self.params.initial_state)
        title = r'${},\ {}$'.format(title1, title2)
        ax.set_title(title)

        ax.set_xlim(left = -0.5, right = len(d) - 0.5)
        ax.set_ylim(bottom = 0, top = 1.1)
        ax.set_xlabel('$j$')
        ax.set_ylabel(ylabel)

        b = ax.bar(x, d, align = 'center')

        fig.savefig(filename)
        plt.close(fig)

    def PlotRho(self, t, mode = 'subspace'):
        """
        Return matrix plot of rho at time t

        """
        title1 = r'\rho'
        rho = self.rho[t]

        if (mode == 'transformed'):
            rho = self.space.ToExcBasis(self.Rho(t))
            title1 = r'\rho^{\rm ex}'

        if (mode == 'full'):
            rho = self.Rho(t)

        title2 = 't = {:7.2f}'.format(float(self.IterationTime(t)))
        title3 = r'\rho(0) = |{0}\rangle \langle {0}|'.format(
                 self.params.initial_state)
        title = r'${},\ {},\ {}$'.format(title1, title2, title3)

        plot_object = matrix_plot(matrix(abs(array(rho))),
            cmap = 'spectral',
            vmin = 0,
            vmax = 1,
            colorbar = True,
            title = title)
        return plot_object

    def PlotSink(self):
        """
        Returns line plot of rho_sink[1,1]

        """
        l = zip(self.TimeList(), self.rho_sink11)
        legend_label = r'$\rho_{1,1}^{\rm sink}(t)$'

        plot_object = line(l,
            ymin = 0,
            ymax = 1,
            color = 'red',
            tick_formatter = 'latex',
            axes_labels = ['$t$', '$P$'],
            legend_label = legend_label)

        plot_object.set_legend_options(back_color = 'white')

        return plot_object

    def ShowSink(self):
        """
        Shows line plot of sink subsystem density matrix elem rho_sink[1,1]

        """
        html('<h2>Sink subsystem density matrix elem</h2>')
        html('Conductivity = {:.4f}'.format(float(self.Conductivity())))
        show(self.PlotSink())
