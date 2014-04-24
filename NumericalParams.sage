class NumericalParams:
    """
    Initializes, stores and prints current values of computation parameters
    """

    def __init__(self, unitary = True):
        self.unitary = unitary

        self.alpha = RDF(1)
        self.beta = RDF(1)
        self.omega_a = RDF(1)
        self.omega_c = RDF(1)
        self.time_steps = 1000
        self.time_end = RDF(40)

        if not unitary:
            self.gamma = RDF(0.1)

    def ShowHTML(self):
        """
        Print values of parameters

        Works inside notebook interface
        """

        l = []

        l.append((r'\alpha', self.alpha))
        l.append((r'\beta', self.beta))

        if not self.unitary:
            l.append((r'\gamma', self.gamma))

        l.append((r'\omega_a', self.omega_a))
        l.append((r'\omega_c', self.omega_c))
        l.append((r'n_t', self.time_steps))
        l.append((r't_{\rm end}', self.time_end)) 

        html_vars(l)
