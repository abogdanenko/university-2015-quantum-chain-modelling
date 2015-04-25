class ConductivityLongChain(object):
    """
    Computes chain conductivity for chains of different lengths

    """
    def __init__(self, max_chain_len = 3):
        self.max_chain_len = max_chain_len
        self.lengths = range(2, self.max_chain_len + 1)

    def ComputeTimeEvolution(self):
        """
        Computes time evolution and conductivity for each chain length

        """
        self.num_list = []

        for n in self.lengths:
            space = Space(chain_len = n)
            num = NumericalComputation(space)
            num.ComputeTimeEvolution()
            self.num_list.append(num)

        self.conductivity = [num.Conductivity() for num in self.num_list]

    def PlotConductivity(self):
        """
        Returns line plot of chain conductivity

        """
        l = zip(self.lengths, self.conductivity)
        legend_label = r'$\langle\rho_{1,1}^{\rm sink}\rangle$'
        labelx = r'$n$'
        labely = r'$C$'

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
