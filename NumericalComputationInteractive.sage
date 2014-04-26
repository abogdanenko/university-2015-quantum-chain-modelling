import matplotlib.pyplot as plt

class NumericalComputationInteractive(NumericalComputation):
    """
    Extends the base class by adding static and interactive plots
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
            axes_labels = ['t', 'probability'],
            legend_label = 'state_number = {}'.format(state))
    
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
        if basis_e:
            m = self.sym.ToEnergyBasis(m)
        title = 'time = {:7.2f}'.format(float(self.IterationTime(t)))
        plot_object = matrix_plot(m,
            cmap = 'gist_heat',
            norm = 'value',
            figsize = 4,
            title = title)
        return plot_object
    
    def ProbabilityBarChart(self, initial_state, t):
        """
        Returns bar chart of state at time t
        """

        title = 'time = {:7.2f}, initial_state = {}'.format(
            float(num.IterationTime(t)), initial_state)

        rho = self.Rho(initial_state, t)
        d = rho.diagonal()
        row = map(abs, d)

        plot_object = bar_chart(row,
            ymin = 0,
            ymax = 1,
            title = title,
            ticks = 1,
            figsize = 4,
            gridlines = True)
        return plot_object
        
    def SaveUnitaryEvolutionMatrixAnimation(self, filename, basis_e = False):
        """
        Saves matrix plots of U as gif animation
        """

        l = [self.PlotUnitaryEvolutionMatrix(t, basis_e)
            for t in range(self.params.time_steps)]
        animation = animate(l)
        animation.gif(savefile = filename, show_path = True)
    
    def SaveProbabilityBarChartAnimation(self, filename, initial_state):
        """
        Saves bar chart of state as gif animation
        """

        l = [self.ProbabilityBarChart(initial_state, t)
            for t in range(self.params.time_steps)]

        animation = animate(l)
        animation.gif(savefile = filename, show_path = True)
    
    def InteractiveParamsSet(self):
        """
        Returns inner function
        
        Prepares default values for inner function
        """

        defaults = NumericalParams(self.sym.unitary)
        alpha_box = input_box(defaults.alpha, label = r'$\alpha = $', type = RDF)
        beta_box = input_box(defaults.beta, label = r'$\beta = $', type = RDF)
        omega_a_box = input_box(defaults.omega_a, label = r'$\omega_a = $', type = RDF)
        omega_c_box = input_box(defaults.omega_c, label = r'$\omega_c = $', type = RDF)
        time_steps_box = input_box(defaults.time_steps, label = '$n_t$', type = Integer)
        time_end_box = input_box(defaults.time_end, label = r'$t_{\rm end}$', type = RDF)

        if self.sym.unitary:
            gamma_box = None
        else:
            gamma_box = input_box(defaults.gamma, label = r'$\gamma = $', type = RDF)
 
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
    
        def inner(initial_state = self.StateSlider(), t = self.TimeSlider()):

            """
            Displays state at time t

            Should be passed to interact()
            """

            html('<h2>State at a given time</h2>')
            show(self.ProbabilityBarChart(initial_state, t))
    
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

            if (basis_e):
                html('<h2>Time evolution matrix in energy basis '
                    'at a given time</h2>')
            else:
                html('<h2>Time evolution matrix at a given time</h2>')
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

            html('<h2>Time evolution of state</h2>')
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

            html('<h2>Time evolution of state vector {}</h2>'.format(initial_state))
            if self.sym.unitary:
                E = energy(initial_state)
                l = self.PlotStates(e_states(E), initial_state)
                show(sum(l))
                show(graphics_array(l, 3, 2))
            else:
                l = self.PlotStates(states_list, initial_state)
                show(sum(l))
                show(graphics_array(l, 8, 2), figsize = [10, 20])

        return inner    

    def EigenVectorsPlot(self):
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
            colorbar = True,
            title = 'Eigen vectors of H (in columns, grouped by E)')
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
            axes_labels = ['t', r'$d(\rho, D)$'],
            legend_label = 'initial_state = {}'.format(initial_state))
    
        return plot_object
    
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
            show(self.PlotDiagDist(initial_state))

        return inner    

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
            axes_labels = ['t', 'probability'],
            legend_label = 'state_number = {}'.format(state))
    
        return plot_object
    
    def PlotReducedStates1(self, initial_state):
        """
        Returns a list of line plots of reduced states
        """

        colors = rainbow(4)
        return [self.PlotReducedState1(s, initial_state, colors[s])
            for s in range(4)]
    
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

            html(('<h2>Time evolution of reduced state vector'
                ' (initial state is {})</h2>').format(initial_state))
            l = self.PlotReducedStates1(initial_state)
            show(sum(l))
            show(graphics_array(l, 2, 2))

        return inner    

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
            ymax = 2,
            axes_labels = ['t', 'entropy'],
            legend_label = 'initial state = {}'.format(initial_state))
    
        return plot_object
    
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
            show(self.PlotEntropy1(initial_state))

        return inner    

    def PlotRho(self, initial_state, t, basis_e = False):
        """
        Return matrix plot of rho at time t
        """

        rho = self.Rho(initial_state, t)
        if basis_e:
            rho = self.sym.ToEnergyBasis(rho)
        title = 'time = {:7.2f}'.format(float(self.IterationTime(t)))
        plot_object = matrix_plot(matrix(abs(array(rho))),
            cmap = 'spectral',
            vmin = 0,
            vmax = 1,
            colorbar = True,
            title = title)
        return plot_object
    
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

            html('<h2>Density matrix at a given time</h2>')
            show(self.PlotRho(initial_state, t, basis_e))
    
        return inner
    
