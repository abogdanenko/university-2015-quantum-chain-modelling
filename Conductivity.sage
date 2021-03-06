from numpy import average

class Conductivity(object):
    """
    Computes chain conductivity, makes sink plots for diff. values of params

    """
    def __init__(self, num):
        self.num = num
        self.param = 'beta'
        # pretty output formatting
        param_list = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 10.0]
        self.param_list = [RDF(x) for x in param_list]

    def InteractiveParamPicker(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        html('<h2>Please pick parameter below</h2>')
        param_picker = selector(
            values = [r'$\beta$', r'$\gamma_s$', r'$\gamma_d$'],
            label = 'Parameter:', buttons = True)

        def inner(param = param_picker):
            """
            Set param to the one chosen by user

            Should be passed to interact()

            """
            html('You have chosen: ' + param)
            # strip dollar signs and backslash
            self.param = param[2:-1]

        return inner

    def ComputeTimeEvolution(self):
        """
        Computes time evolution for each parameter value

        """
        self.rho_sink11_list = []
        self.conductivity = []
        for param in self.param_list:
            if self.param == 'beta':
                self.num.params.beta = param
            elif self.param == 'gamma_s':
                self.num.params.gamma_s = param
            else:
                self.num.params.gamma_d = param
            self.num.ParamsChanged()
            self.num.ComputeTimeEvolution()
            l = self.num.rho_sink11
            self.rho_sink11_list.append(l)
            c = self.num.Conductivity()
            self.conductivity.append(c)

    def PlotSink(self, param_index, color = 'red'):
        """
        Returns line plot of sink subsystem matrix elem

        """
        l = zip(self.num.TimeList(), self.rho_sink11_list[param_index])

        param = '{}'.format(self.param_list[param_index])
        label1 = r'\rho_{1,1}^{\rm sink}(t)'
        label2 = r'\{} = {}'.format(self.param, param)
        legend_label = r'${},\ {}$'.format(label1, label2)
        plot_object = line(l,
            ymin = 0,
            ymax = 1,
            color = color,
            tick_formatter = 'latex',
            axes_labels = ['$t$', '$P$'],
            legend_label = legend_label)

        plot_object.set_legend_options(back_color = 'white')

        return plot_object

    def PlotSinks(self):
        """
        Returns multi-line plot of rho_sink[1,1] with different param values

        """
        n = len(self.param_list)
        colors = rainbow(n)

        plot = sage.plot.graphics.Graphics()
        for i in range(n):
            plot += self.PlotSink(i, colors[i])

        return plot

    def ShowSink(self):
        """
        Shows multi-line plot of rho_sink[1,1] with different param values

        """
        html('<h2>Sink subsystem density matrix elem</h2>')
        show(self.PlotSinks())

    def InteractiveSink(self):
        """
        Returns inner function

        Prepares default values for inner function

        """
        html('<h2>Sink subsystem density matrix elem</h2>')

        param_value_picker = selector(values = self.param_list,
            label = '{}:'.format(self.param))

        def inner(param = param_value_picker):
            """
            Shows line plot of sink subsystem matrix elem
            User specifies param

            Should be passed to interact()

            """
            index = self.param_list.index(param)
            html('Conductivity = {:.4f}'.format(float(
                self.conductivity[index])))
            show(self.PlotSink(index))

        return inner

    def PlotConductivity(self):
        """
        Returns line plot of chain conductivity

        """
        l = zip(self.param_list, self.conductivity)
        legend_label = r'$\langle\rho_{1,1}^{\rm sink}\rangle$'
        labelx = r'$\{}$'.format(self.param)
        labely = '$C$'

        plot_object = line(l,
            ymin = 0,
            ymax = 1,
            xmin = 0,
            color = 'red',
            tick_formatter = 'latex',
            axes_labels = [labelx, labely],
            legend_label = legend_label)

        plot_object.set_legend_options(back_color = 'white')

        return plot_object

    def ShowConductivity(self):
        """
        Shows line plot of chain conductivity

        """
        html('<h2>Chain conductivity</h2>')
        show(self.PlotConductivity())
