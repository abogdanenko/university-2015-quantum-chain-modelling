def exc_number(i):
    """
    Returns number of ones in binary notation of i

    Which is coincidentally total number of excitations
    """
    return Integer(i).bits().count(1)

def e_states(E):
    """
    Returns a list of states with total number of excitations E

    """
    return [i for i in states_list if exc_number(i) == E]

def vec2dm(psi):
    """
    Returns density matrix corresponding to state vector psi

    """
    c = psi.column()
    return c * c.conjugate_transpose()

def full_basis_state(state):
    """
    Returns basis state number ``state``

    """
    states_count = 2 ** 5
    psi = vector(CDF, states_count)
    psi[state * 2] = 1
    return psi

def coordinate_change_matrix():
    """
    Returns coordinate change matrix to excitation basis

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

def partial_trace1(rho):
    """
    Returns reduced density matrix of subsystem 1

    """
    rho1 = matrix(CDF, 4)
    for i1i2 in range(4):
        for j1j2 in range(4):
             s = 0
             for i3i4 in range(4):
                 s += rho[i1i2 * 4 + i3i4, j1j2 * 4 + i3i4]

             rho1[i1i2, j1j2] = s
    return rho1

def partial_trace_chain(rho_full):
    """
    Returns reduced density matrix of chain subsystem

    """
    rho = matrix(CDF, 16)
    for i in range(16):
        for j in range(16):
             s = 0
             for k in range(2):
                 s += rho_full[i * 2 + k, j * 2 + k]

             rho[i, j] = s
    return rho

def partial_trace_sink(rho_full):
    """
    Returns reduced density matrix of sink subsystem

    """
    rho = matrix(CDF, 2)
    for i in range(2):
        for j in range(2):
             s = 0
             for k in range(16):
                 s += rho_full[i + 2 * k, j + 2 * k]

             rho[i, j] = s
    return rho
