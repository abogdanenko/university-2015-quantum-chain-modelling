import matplotlib.pyplot as plt

class SymbolicComputationPlots(SymbolicComputationBase):
    """
    Makes plots

    """
    def EigenValuesBarChart(self, alpha, beta):
        """
        Displays bar chart of eigenvalues grouped by total number of excitations

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
        fig.savefig('EigenValuesBarChartPyPlt')
