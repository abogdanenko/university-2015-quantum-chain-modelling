import matplotlib.pyplot as plt

class NumericalComputationPlots(NumericalComputationBase):
    """
    Makes plots

    """
    def PlotState(self, state, color = 'red'):
        """
        Returns line plot of the state

        """
        l = []
        for t in range(self.params.time_steps):
            x = self.IterationTime(t)

            rho = self.Rho(t)
            y = abs(rho[state, state])

            point = (x, y)
            l.append(point)


        plot_object = line(l,
            ymin = 0,
            ymax = 1,
            color = color,
            tick_formatter = 'latex',
            axes_labels = ['$t$', '$P$'],
            legend_label = r'$\rho_{{{0},{0}}}(t)$'.format(state))

        plot_object.set_legend_options(back_color = 'white')

        return plot_object

    def PlotStates(self, states):
        """
        Returns a list of line plots of the states

        """
        colors = rainbow(len(states))

        return [self.PlotState(s, c) for s, c \
            in zip(states, colors)]

    def ShowStates(self):
        """
        Shows multi-line plot of state evolution

        """
        html('<h2>Density matrix diagonal elems</h2>')
        l = self.PlotStates(states_list)
        show(sum(l))

    def ProbabilityBarChart(
            self,
            t,
            basis_e = False,
            filename = 'ProbabilityBarChart'):
        """
        Returns bar chart of state at time t

        """
        rho = self.Rho(t)
        if basis_e:
            rho = self.sym.ToExcBasis(rho)
        d = map(abs, rho.diagonal())

        colors = []
        if basis_e:
            ylabel = r'$\rho_{j,j}^{\rm ex}$'
            for E in exc_list:
                colors.extend(block_sizes[E] * [exc_number_rainbow[E]])
        else:
            ylabel = r'$\rho_{j,j}$'
            for state in states_list:
                colors.append(exc_number_rainbow[exc_number(state)])

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        title1 = 't = {:7.2f}'.format(float(num.IterationTime(t)))
        title2 = r'\rho(0) = |{0}\rangle \langle {0}|'.format(
                 self.params.initial_state)
        title = r'${},\ {}$'.format(title1, title2)
        ax.set_title(title)

        ax.set_xlim(left = -0.5, right = states_count - 0.5)
        ax.set_ylim(bottom = 0, top = 1.1)
        ax.set_xlabel('$j$')
        ax.set_ylabel(ylabel)

        ax.set_xticks(states_list)
        ax.set_xticklabels(states_list)

        b = ax.bar(states_list, d, color = colors, align = 'center')

        indices = left_indices if basis_e else first_indices
        labels = [r'$N_{{\rm ex}} = {}$'.format(E) for E in exc_list]
        points = [b[indices[E]] for E in exc_list]
        ax.legend(points, labels, loc = 'upper center')
        fig.savefig(filename)
        plt.close(fig)

    def PlotRho(self, t, basis_e = False):
        """
        Return matrix plot of rho at time t

        """
        rho = self.Rho(t)
        title1 = r'\rho'
        if basis_e:
            rho = self.sym.ToExcBasis(rho)
            title1 = r'\rho^{\rm ex}'

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

    def PlotSink(self, i, j):
        """
        Returns line plot of sink subsystem matrix elem

        """
        l = []
        for t in range(self.params.time_steps):
            x = self.IterationTime(t)

            rho = self.rho_sink_list[t]

            y = abs(rho[i, j])

            point = (x, y)
            l.append(point)

        index = r'{},{}'.format(i, j)
        legend_label = r'$\rho_{' + index + r'}^{\rm sink}(t)$'
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
        show(self.PlotSink(1, 1))
