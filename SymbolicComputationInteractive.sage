class SymbolicComputationInteractive(SymbolicComputationPlots):
    """
    Makes interactive plots

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
            html(r'<h3>Eigen values of H (grouped by $N_{\rm ex}$)</h3>')
            if name == r'$\alpha$':
                p = self.EigenValuesLinePlot(
                    omega = defaults.omega_a,
                    alpha = alphabeta)
            else:
                p = self.EigenValuesLinePlot(
                    omega = defaults.omega_a,
                    beta = alphabeta)

            show(p)

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
            html(r'<h3>Eigen values of H (grouped by $N_{\rm ex}$)</h3>')
            self.EigenValuesBarChart(alpha, beta)

        return inner
