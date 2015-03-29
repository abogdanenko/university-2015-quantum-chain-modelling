class Conductivity(object):
    """
    Computes chain conductivity, makes sink plots for diff. values of params

    """
    def __init__(self, num):
        self.num = num
        # pretty output formatting
        beta_list = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 10.0]
        self.beta_list = [RDF(x) for x in beta_list]

    def GetRhoSink11(self):
        """
        Return list of rho_sink[1,1] for each moment of time

        """
        l = []
        for t in range(self.num.params.time_steps):
            rho = self.num.rho_sink_list[t]
            y = abs(rho[1, 1])
            l.append(y)
        return l

    def ComputeTimeEvolution(self):
        """
        Computes time evolution for each beta

        """
        self.rho_sink_11_list = []
        for beta in self.beta_list:
            self.num.params.beta = beta
            self.num.InitOperators()
            self.num.ComputeTimeEvolution()
            l = self.GetRhoSink11()
            self.rho_sink_11_list.append(l)

    def PlotSink(self, beta_index, color = 'red'):
        """
        Returns line plot of sink subsystem matrix elem

        """
        l = []
        for t in range(self.num.params.time_steps):
            x = self.num.IterationTime(t)
            y = self.rho_sink_11_list[beta_index][t]
            point = (x, y)
            l.append(point)

        beta = '{}'.format(self.beta_list[beta_index])
        legend_label = r'$\rho_{1,1}^{\rm sink}(t),\ \beta = ' + beta + '$'
        plot_object = line(l,
            ymin = 0,
            ymax = 1,
            color = color,
            tick_formatter = 'latex',
            axes_labels = ['$t$', '$P$'],
            legend_label = legend_label)

        plot_object.set_legend_options(back_color = 'white')

        return plot_object

    def ShowSink(self):
        """
        Shows multi-line plot of rho_sink[1,1] with different beta

        """
        html('<h2>Sink subsystem density matrix elem</h2>')
        n = len(self.beta_list)
        colors = rainbow(n)

        plot = sage.plot.graphics.Graphics()
        for i in range(n):
            plot += self.PlotSink(i, colors[i])

        return plot

    def InteractiveSink(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        html('<h2>Sink subsystem density matrix elem</h2>')
        def inner(beta = self.beta_list):
            """
            Shows line plot of sink subsystem matrix elem
            User specifies beta

            Should be passed to interact()

            """
            index = self.beta_list.index(beta)
            show(self.PlotSink(index))

        return inner
