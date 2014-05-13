class NumericalComputationInteractive(object):
    """
    Makes interactive plots for class NumericalComputation

    """
    def InteractiveParamsSet(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        defaults = NumericalParams(self.sym.unitary)

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

        if self.sym.unitary:
            gamma_box = None
        else:
            gamma_box = input_box(
                defaults.gamma,
                label = r'$\gamma = $',
                type = RDF)

        # todo: refactor two function definitions into one

        def inner(
                alpha = alpha_box,
                beta = beta_box,
                omega_a = omega_a_box,
                omega_c = omega_c_box,
                time_steps = time_steps_box,
                time_end = time_end_box,
                auto_update = False):
            """
            Sets parameters interactively, display new parameters in html

            Should be passed to interact()

            """
            self.params.alpha = alpha
            self.params.beta = beta
            self.params.omega_a = omega_a
            self.params.omega_c = omega_c
            self.params.time_end = time_end
            self.params.time_steps = time_steps
            self.InitOperators()
            self.ComputeEigenVectors()
            self.params.ShowHTML()

        def inner_not_unitary(
                alpha = alpha_box,
                beta = beta_box,
                gamma = gamma_box,
                omega_a = omega_a_box,
                omega_c = omega_c_box,
                time_steps = time_steps_box,
                time_end = time_end_box,
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
            self.InitOperators()
            self.ComputeEigenVectors()
            self.params.ShowHTML()

        if self.sym.unitary:
            return inner
        else:
            return inner_not_unitary

    def StateSlider(self):
        return slider(states_list, default = 1)

    def TimeSlider(self):
        return slider(range(self.params.time_steps))

    def InteractiveProbabilityBarChart(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(
                initial_state = self.StateSlider(),
                t = self.TimeSlider(),
                basis_e = [True, False]):
            """
            Displays state at time t

            Should be passed to interact()

            """
            self.ProbabilityBarChart(initial_state, t, basis_e)

        return inner

    def InteractiveUnitaryEvolutionMatrix(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(
                t = self.TimeSlider(),
                basis_e = [True, False]):
            """
            Displays time evolution matrix at time t

            Should be passed to interact()

            """
            html('<h2>Time evolution matrix</h2>')
            show(self.PlotUnitaryEvolutionMatrix(t, basis_e))

        return inner

    def InteractiveState(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(
                initial_state = self.StateSlider(),
                state = self.StateSlider(),
                auto_update = False):
            """
            Shows line plot of state vector component
            User specifies initial state and component to plot

            Should be passed to interact()

            """
            html(r'$\rho(0) = |{0}\rangle \langle {0}|$'.format(
                initial_state))
            show(self.PlotState(state, initial_state))

        return inner

    def InteractiveStates(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(initial_state = self.StateSlider(), auto_update = False):
            """
            Shows multi-line plot and then multiple line plots of state evol.

            User specifies initial state

            """
            html(r'$\rho(0) = |{0}\rangle \langle {0}|$'.format(
                initial_state))
            if self.sym.unitary:
                E = exc_number(initial_state)
                l = self.PlotStates(e_states(E), initial_state)
                show(sum(l))
                show(graphics_array(l, 3, 2))
            else:
                l = self.PlotStates(states_list, initial_state)
                show(sum(l))
                show(graphics_array(l, 8, 2), figsize = [10, 20])

        return inner

    def InteractiveDiagDist(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(initial_state = self.StateSlider()):
            """
            Shows line plot of DiagDist(initial_state, t)

            User specifies initial state

            """
            html('<h2>Distance to diagonal matrices</h2>')
            html(r'$\rho(0) = |{0}\rangle \langle {0}|$'.format(
                initial_state))
            show(self.PlotDiagDist(initial_state))

        return inner

    def InteractiveReducedStates1(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(initial_state = self.StateSlider()):
            """
            Shows multi-line plot and then multiple line plots of state evol.

            User specifies initial state

            """
            html(r'$\rho(0) = |{0}\rangle \langle {0}|$'.format(
                initial_state))
            l = self.PlotReducedStates1(initial_state)
            show(sum(l))
            show(graphics_array(l, 2, 2))

        return inner

    def InteractiveEntropy1(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(initial_state = self.StateSlider()):
            """
            Shows line plot of entropy

            User specifies initial state

            """
            html('<h2>Entropy of subsystem 1</h2>')
            html(r'$\rho(0) = |{0}\rangle \langle {0}|$'.format(
                initial_state))
            show(self.PlotEntropy1(initial_state))

        return inner

    def InteractiveRho(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        def inner(
                initial_state = self.StateSlider(),
                t = self.TimeSlider(),
                basis_e = [True, False]):
            """
            Displays density matrix at time t

            Should be passed to interact()

            """
            show(self.PlotRho(initial_state, t, basis_e))

        return inner

    def ShowEigenVectorsPlot(self):
        html(r'<h3>Eigen vectors of H (in columns, grouped by $N_{\rm ex}$)</h3>')
        show(self.PlotEigenVectors())
