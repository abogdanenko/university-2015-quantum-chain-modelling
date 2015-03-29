class Conductivity(object):
    """
    Computes chain conductivity, makes sink plots for diff. values of params

    """
    def __init__(self, num):
        self.num = num
        self.beta_list = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 10.0]

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

    def ComputeTimeEvolutionBeta(self):
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
