class NumericalComputationInteractive(NumericalComputationPlots):
    """
    Makes interactive plots

    """
    def InteractiveParamsSet(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        html('<h2>Set parameters, number of points, time interval '
            'and initial condition</h2>')

        defaults = NumericalParams()

        alpha_box = input_box(
            defaults.alpha,
            label = r'$\alpha = $',
            type = RDF)

        beta_box = input_box(
            defaults.beta,
            label = r'$\beta = $',
            type = RDF)

        omega_a_box = input_box(
            defaults.omega_a,
            label = r'$\omega_a = $',
            type = RDF)

        omega_c_box = input_box(
            defaults.omega_c,
            label = r'$\omega_c = $',
            type = RDF)

        time_steps_box = input_box(
            defaults.time_steps,
            label = '$n_t$',
            type = Integer)

        time_end_box = input_box(
            defaults.time_end,
            label = r'$t_{\rm end}$',
            type = RDF)

        initial_state_slider = self.space.StateSlider(
            default = defaults.initial_state)

        gamma_box = input_box(
            defaults.gamma,
            label = r'$\gamma = $',
            type = RDF)

        def inner(
                alpha = alpha_box,
                beta = beta_box,
                gamma = gamma_box,
                omega_a = omega_a_box,
                omega_c = omega_c_box,
                time_steps = time_steps_box,
                time_end = time_end_box,
                initial_state = initial_state_slider,
                auto_update = False):
            """
            Sets parameters interactively, display new parameters in html

            Should be passed to interact()

            """
            self.params.alpha = alpha
            self.params.beta = beta
            self.params.gamma = gamma
            self.params.omega_a = omega_a
            self.params.omega_c = omega_c
            self.params.time_end = time_end
            self.params.time_steps = time_steps
            self.params.initial_state = initial_state
            self.ParamsChanged()
            self.params.ShowHTML()

        return inner

    def TimeSlider(self):
        """
        Returns slider to select point in time

        """
        return slider(range(self.params.time_steps))

    def InteractiveProbabilityBarChart(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        html('<h2>Density matrix diagonal</h2>')
        def inner(
                t = self.TimeSlider(),
                mode = ['subspace', 'transformed', 'full']):
            """
            Displays state at time t

            Should be passed to interact()

            """
            self.ProbabilityBarChart(t, mode)

        return inner

    def InteractiveState(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        html('<h2>Density matrix diagonal elem</h2>')
        def inner(
                state = self.space.StateSlider(),
                auto_update = False):
            """
            Shows line plot of state vector component
            User specifies initial state and component to plot

            Should be passed to interact()

            """
            show(self.PlotState(state))

        return inner

    def InteractiveRho(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        html('<h2>Density matrix</h2>')
        def inner(
                t = self.TimeSlider(),
                mode = ['subspace', 'transformed', 'full']):
            """
            Displays density matrix at time t

            Should be passed to interact()

            """
            show(self.PlotRho(t, mode))

        return inner
