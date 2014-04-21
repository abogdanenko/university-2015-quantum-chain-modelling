class SymbolicComputationInteractive(SymbolicComputation):
    """
    Extends the base class by adding static and interactive plots
    """

    def InteractiveEigenValuesLinePlot(self):
        """
        Returns inner function
        
        Prepares default values for inner function
        """

        values = [RealField(10)(k / 10.0) for k in range(41)]
        defaults = NumericalParams()

        alphabeta_slider = slider(
            values,
            default = defaults.alpha,
            label = 'Value of the fixed variable: ')

        name_picker = ('Fix one variable:', [r'$\alpha$', r'$\beta$'])

        def inner(name = name_picker, alphabeta = alphabeta_slider):
            """
            Shows line plot of eigen values
            User specifies a variable (alpha or beta) and a value

            Should be passed to interact()
            """

            def eigenvalues(E):
                kwds = {'omega': defaults.omega_a}
                key = 'alpha' if name == r'$\alpha$' else 'beta'
                kwds[key] = alphabeta
                return [value.subs(**kwds)
                    for value in self.eigenvalues_blocks[E]]

            xlabel = r'$\alpha$' if name == r'$\beta$' else r'$\beta$'

            plots = [plot(
                    eigenvalues(E),
                    ymin = -6,
                    ymax = 8,
                    xmin = 0,
                    xmax = 4,
                    color = energy_rainbow[E],
                    legend_label = 'E = {}'.format(E),
                    axes_labels = [xlabel, 'Energy'],
                    title = 'Eigen values of H (grouped by E)')
                for E in energy_list]
            show(sum(plots))

        return inner

    def EigenValuesBarChart(self, alpha, beta):
        """
        Displays bar chart of eigenvalues grouped by energy
        """
        
        values = []
        colors = []
        for E in energy_list:
            for value in self.eigenvalues_blocks[E]:
                v = value.subs(
                    omega = NumericalParams().omega_a,
                    alpha = alpha,
                    beta = beta)
                values.append(v)
                colors.append(energy_rainbow[E])


        fig = plt.figure()
        ax = fig.add_subplot(111r)
        ax.set_title('Eigen values of H (grouped by E)')
        ax.set_xticklabels([])
        ax.set_xlim(left = 0, right = states_count)
        ax.set_ylabel('Energy')
        b = ax.bar(states_list, values, color = colors, linewidth = 4)
        points = [b[left_indices[E]] for E in energy_list]
        labels = ['E = {}'.format(E) for E in energy_list]
        ax.legend(points, labels, loc = 'upper left')
        ax.grid(True)
        fig.savefig('EigenValuesBarChartPyPlt')        

    def InteractiveEigenValuesBarChart(self):
        """
        Returns inner function
        
        Prepares default values for inner function
        """

        values = [RealField(10)(k / 10.0) for k in range(41)]
        defaults = NumericalParams()

        alpha_slider = slider(
            values,
            default = defaults.alpha,
            label = r'$\alpha: $')

        beta_slider = slider(
            values,
            default = defaults.beta,
            label = r'$\beta: $')

        def inner(alpha = alpha_slider, beta = beta_slider):
            """
            Shows bar chart of eigen values
            User specifies alpha and beta

            Should be passed to interact()
            """

            self.EigenValuesBarChart(alpha, beta)

        return inner
