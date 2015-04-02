class NumericalParams(object):
    """
    Initializes, stores and prints current values of computation parameters

    """
    def __init__(self):
        self.alpha = RDF(1)
        self.beta = RDF(1)
        self.omega_a = RDF(1)
        self.omega_c = RDF(1)
        self.time_steps = 500
        self.time_end = RDF(40)
        self.gamma = RDF(0.5)
        self.initial_state = 8

    def Dt(self):
        return RDF(self.time_end / self.time_steps)

    def ShowHTML(self):
        """
        Print values of parameters

        Works inside notebook interface

        """
        l = []

        l.append((r'\alpha', self.alpha))
        l.append((r'\beta', self.beta))
        l.append((r'\gamma', self.gamma))
        l.append((r'\omega_a', self.omega_a))
        l.append((r'\omega_c', self.omega_c))
        l.append((r'n_t', self.time_steps))
        l.append((r't_{\rm end}', self.time_end))

        html_vars(l)
        html(r'$\rho(0) = |{0}\rangle \langle {0}|$'.format(self.initial_state))
