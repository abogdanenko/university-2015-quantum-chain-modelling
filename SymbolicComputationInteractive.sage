class SymbolicComputationInteractive(object):
    """
    Makes interactive plots for class SymbolicComputation

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
                    color = exc_number_rainbow[E],
                    legend_label = r'$N_{{\rm ex}}$ = {}'.format(E),
                    axes_labels = [xlabel, '$E$'],
                    title = r'Eigen values of H (grouped by $N_{\rm ex}$)')
                for E in exc_list]
            show(sum(plots))

        return inner

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
