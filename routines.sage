def energy(i):
    """
    Returns number of ones in binary notation of i

    Which is coincidentally total number of excitations
    """
    return Integer(i).bits().count(1)

def e_states(E):
    """
    Returns a list of states with total number of excitations E

    """
    return [i for i in states_list if energy(i) == E]

def vec2dm(psi):
    """
    Returns density matrix corresponding to state vector psi

    """
    c = psi.column()
    return c * c.conjugate_transpose()

def basis_state(state):
    """
    Returns basis state number ``state``

    """
    psi = vector(CDF, states_count)
    psi[state] = 1
    return psi

def coordinate_change_matrix():
    """
    Returns coordinate change matrix to energy basis

    """
    T = matrix(states_count)
    rows = []
    for E in exc_list:
        rows.extend(e_states(E))
    for j in states_list:
        i = rows[j]
        T[i, j] = 1
    return T

def bits(n, min_width = 1):
    """
    Returns binary representation of n as a list of ones and zeroes

    The returned list starts with a major binary digit

    """
    l = Integer(n).bits()
    l.reverse()
    padded_list = [0] * (min_width - len(l)) + l
    return padded_list

def html_vars(l):
    """
    Prints variables and values from a list of pairs

    Works with notebook interface

    """
    for x in l:
        html('${} = {}$'.format(x[0], latex(x[1])))

def flatten_eigen(ev):
    """
    Returnes a list of eigenvalue, eigenvector pairs

    ev should be in fmt. returned by Matrix_symbolic_dense.eigenvectors_right()

    """
    l = []
    for value, vectors, multiplicity in ev:
        for vector in vectors:
            t = (value, vector)
            l.append(t)
    return l

def ev_flat_sorted(A):
    """
    Returnes a sorted list of eigenvalue, eigenvector pairs of matrix A

    """
    ev = A.eigenvectors_right()
    ev = flatten_eigen(ev)
    return sorted(ev)

def show_eigen_html(ev):
    """
    Prints eigenvalues and eigenvectors

    Works in notebook interface.

    """
    def sparse_format(x):
        epsilon = 1e-10
        if (abs(x) < epsilon):
            return ' ' * 7
        else:
            return '{:7.2f}'.format(float(x))

    if (ev[0][0].base_ring() == SR):
        html.table(ev, header = ['eigenvalue', 'eigenvector'])
    else:
        l = []
        for value, vector in ev:
            string_list = [sparse_format(x) for x in vector]
            l.append(('{:7.2f}'.format(float(value)),
                '({})'.format(''.join(string_list))))
        html.table(l, header = ['eigenvalue', 'eigenvector'])

def partial_trace1(rho):
    """
    Returnes reduced density matrix of subsystem 1

    """
    rho1 = matrix(CDF, 4)
    for i1i2 in range(4):
        for j1j2 in range(4):
             s = 0
             for i3i4 in range(4):
                 s += rho[i1i2 * 4 + i3i4, j1j2 * 4 + i3i4]

             rho1[i1i2, j1j2] = s
    return rho1

def xlnx(x):
    """
    Returnes x * ln(x) or 0 if x is near 0

    """
    epsilon = 1e-14
    if (abs(x) < epsilon):
        return 0
    else:
        return x * ln(x)
