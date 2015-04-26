class NumericalParams(object):
    """
    Initializes, stores and prints current values of computation parameters

    """
    def __init__(self, space):
        self.alpha = RDF(0.7)
        self.beta = RDF(1.3)
        self.omega_a = RDF(0.8)
        self.omega_c = RDF(1.2)
        self.gamma_s = RDF(3.0)
        self.gamma_d = RDF(0.5)
        self.time_steps = 500
        self.time_end = RDF(40)
        # second to last state from subspace 1 has exciton 1 in first cavity
        self.initial_state = space.subspaces[1].states[-2].index

    def __repr__(self):
        s = ''
        s += 'alpha = {}; '.format(self.alpha)
        s += 'beta = {}; '.format(self.beta)
        s += 'omega_a = {}; '.format(self.omega_a)
        s += 'omega_c = {}; '.format(self.omega_c)
        s += 'gamma_s = {}; '.format(self.gamma_s)
        s += 'gamma_d = {}; '.format(self.gamma_d)
        s += 'time_steps = {}; '.format(self.time_steps)
        s += 'time_end = {}; '.format(self.time_end)
        s += 'initial_state = {}'.format(self.initial_state)
        return s

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
        l.append((r'\gamma_s', self.gamma_s))
        l.append((r'\gamma_d', self.gamma_d))
        l.append((r'\omega_a', self.omega_a))
        l.append((r'\omega_c', self.omega_c))
        l.append((r'n_t', self.time_steps))
        l.append((r't_{\rm end}', self.time_end))

        html_vars(l)
        html(r'$\rho(0) = |{0}\rangle \langle {0}|$'.format(self.initial_state))
