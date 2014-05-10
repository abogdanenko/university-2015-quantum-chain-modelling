import matplotlib.pyplot as plt

class SymbolicComputationPlots(object):
    """
    Makes plots for class SymbolicComputation
    """

    def EigenValuesBarChart(self, alpha, beta):
        """
        Displays bar chart of eigenvalues grouped by energy
        """
        
        values = []
        colors = []
        for E in energy_list:
            for value in self.eigenvalues_blocks[E]:
                v = value.subs(
                    omega = NumericalParams().omega_a,
                    alpha = alpha,
                    beta = beta)
                values.append(v)
                colors.append(energy_rainbow[E])


        fig = plt.figure()
        ax = fig.add_subplot(111r)
        ax.set_title('Eigen values of H (grouped by E)')
        ax.set_xticklabels([])
        ax.set_xlim(left = 0, right = states_count)
        ax.set_ylabel('Energy')
        b = ax.bar(states_list, values, color = colors, linewidth = 4)
        points = [b[left_indices[E]] for E in energy_list]
        labels = ['E = {}'.format(E) for E in energy_list]
        ax.legend(points, labels, loc = 'upper left')
        ax.grid(True)
        fig.savefig('EigenValuesBarChartPyPlt')        
