import matplotlib.pyplot as plt

class SymbolicComputationPlots(SymbolicComputationBase):
    """
    Makes plots

    """
    def EigenValuesLinePlot(self, **kwds):
        """
        Returns line plot of eigen values

        kwds must include: (alpha or beta) and omega

        """
        xlabel = r'$\alpha$' if 'beta' in kwds else r'$\beta$'

        plots = []
        for E in exc_list:
            k = range(block_sizes[E])
            k_str = ','.join(map(str, k))
            label = r'$E_{{{}}}^{}$'.format(k_str, E)

            ev_block = self.eigenvalues_blocks[E]
            eigenvalues = [value.subs(**kwds) for value in ev_block]

            p = plot(
                    eigenvalues,
                    ymin = -6,
                    ymax = 8,
                    xmin = 0,
                    xmax = 4,
                    tick_formatter = 'latex',
                    color = exc_number_rainbow[E],
                    legend_label = label,
                    axes_labels = [xlabel, '$E$'])
            p.set_legend_options(back_color = 'white')
            p.set_legend_options(loc = 'upper center')
            p.set_legend_options(ncol = len(exc_list))
            plots.append(p)
        return sum(plots)

    def EigenValuesBarChart(
            self,
            alpha,
            beta,
            filename = 'EigenValuesBarChartPyPlt'):
        """
        Displays bar chart of eigenvalues grouped by total number of excitations

        Use filename parameter to export the plot to file
        """
        values = []
        colors = []
        for E in exc_list:
            for value in self.eigenvalues_blocks[E]:
                v = value.subs(
                    omega = NumericalParams().omega_a,
                    alpha = alpha,
                    beta = beta)
                values.append(v)
                colors.append(exc_number_rainbow[E])

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xticks(states_list)
        ax.set_xlim(left = -0.5, right = states_count - 0.5)
        ax.set_xlabel(r'$(N_{\rm ex}, k)$')

        labels = []
        for E in exc_list:
            for k in range(block_sizes[E]):
                label = r'$({}, {})$'.format(E, k)
                labels.append(label)

        ax.set_xticklabels(labels)
        ax.set_ylabel('$E$')
        b = ax.bar(
            states_list,
            values,
            color = colors,
            linewidth = 4,
            align = 'center')
        points = [b[left_indices[E]] for E in exc_list]
        labels = [r'$E_k^{}$'.format(E) for E in exc_list]
        ax.legend(points, labels, loc = 'upper left')
        ax.grid(True)
        fig.savefig(filename)
