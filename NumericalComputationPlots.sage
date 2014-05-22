import matplotlib.pyplot as plt

class NumericalComputationPlots(NumericalComputationBase):
    """
    Makes plots

    """
    def PlotState(self, state, initial_state = None, color = 'red'):
        """
        Returns line plot of the state

        """
        if (initial_state == None):
            initial_state = state

        l = []
        for t in range(self.params.time_steps):
            x = self.IterationTime(t)

            rho = self.Rho(initial_state, t)
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

    def PlotStates(self, states, initial_state = None):
        """
        Returns a list of line plots of the states

        """
        if (initial_state == None):
            initial_state = states[0]
        colors = rainbow(len(states))

        return [self.PlotState(s, initial_state, c) for s, c \
            in zip(states, colors)]

    def PlotUnitaryEvolutionMatrix(self, t, basis_e = False):
        """
        Return matrix plot of U at time t

        """
        m = self.U[t].apply_map(norm)
        title1 = 'U'
        if basis_e:
            m = self.sym.ToExcBasis(m)
            title1 = r'U^{\rm ex}'

        title2 = 't = {:7.2f}'.format(float(self.IterationTime(t)))
        title = r'${},\ {}$'.format(title1, title2)

        plot_object = matrix_plot(m,
            cmap = 'gist_heat',
            vmin = 0,
            vmax = 1,
            tick_formatter = 'latex',
            colorbar = True,
            title = title)
        return plot_object

    def ProbabilityBarChart(
            self,
            initial_state,
            t,
            basis_e = False,
            filename = 'ProbabilityBarChart'):
        """
        Returns bar chart of state at time t

        """
        rho = self.Rho(initial_state, t)
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
                 initial_state)
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

    def PlotEigenVectors(self):
        """
        Return eigen matrix plot

        todo: make matrix m in NumericalComputation.ComputeEigenVectors

        """
        m = matrix(RDF, states_count)
        column = 0
        row = 0
        for ev in self.ev_blocks:
            for value, vector in ev:
                for i in range(len(vector)):
                    m[row + i, column] = vector[i]
                column += 1
            row += len(vector)

        plot_object = matrix_plot(m,
            cmap = 'bwr',
            vmin = -1,
            vmax = 1,
            tick_formatter = 'latex',
            colorbar = True)
        return plot_object

    def PlotDiagDist(self, initial_state):
        """
        Returns line plot of DiagDist(initial_state, t)

        """
        l = []
        for t in range(self.params.time_steps):
            x = self.IterationTime(t)
            y = self.DiagDist(initial_state, t)
            point = (x, y)
            l.append(point)

        plot_object = line(l,
            ymin = 0,
            ymax = 1,
            tick_formatter = 'latex',
            axes_labels = ['$t$', r'$d(\rho, D)$'])

        return plot_object

    def PlotReducedState1(self, state, initial_state, color = 'red'):
        """
        Returns line plot of reduced state

        """
        l = []
        for t in range(self.params.time_steps):
            x = self.IterationTime(t)
            rho = self.Rho(initial_state, t)
            rho1 = partial_trace1(rho)
            y = abs(rho1[state, state])
            point = (x, y)
            l.append(point)

        plot_object = line(l,
            ymin = 0,
            ymax = 1,
            color = color,
            tick_formatter = 'latex',
            axes_labels = ['$t$', '$P$'],
            legend_label = r'$\rho^1_{{{0},{0}}}(t)$'.format(state))

        plot_object.set_legend_options(back_color = 'white')

        return plot_object

    def PlotReducedStates1(self, initial_state):
        """
        Returns a list of line plots of reduced states

        """
        colors = rainbow(4)
        return [self.PlotReducedState1(s, initial_state, colors[s])
            for s in range(4)]

    def PlotEntropy1(self, initial_state):
        """
        Returns line plot of entropy

        """
        l = []
        for t in range(self.params.time_steps):
            x = self.IterationTime(t)
            y = self.Entropy1(initial_state, t)
            point = (x, y)
            l.append(point)

        plot_object = line(l,
            ymin = 0,
            tick_formatter = 'latex',
            axes_labels = ['$t$', r'$S(\rho^1)$'])

        return plot_object

    def PlotRho(self, initial_state, t, basis_e = False):
        """
        Return matrix plot of rho at time t

        """
        rho = self.Rho(initial_state, t)
        title1 = r'\rho'
        if basis_e:
            rho = self.sym.ToExcBasis(rho)
            title1 = r'\rho^{\rm ex}'

        title2 = 't = {:7.2f}'.format(float(self.IterationTime(t)))
        title3 = r'\rho(0) = |{0}\rangle \langle {0}|'.format(
                 initial_state)
        title = r'${},\ {},\ {}$'.format(title1, title2, title3)

        plot_object = matrix_plot(matrix(abs(array(rho))),
            cmap = 'spectral',
            vmin = 0,
            vmax = 1,
            tick_formatter = 'latex',
            colorbar = True,
            title = title)
        return plot_object
