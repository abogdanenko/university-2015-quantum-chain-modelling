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

def exc_index(i):
    """
    Returns index of state i within its subspace

    """
    E = exc_number(i)
    subspace = e_states(E)
    return subspace.index(i)

def exc_size(i):
    """
    Returns size of subspace to which i belongs

    """
    E = exc_number(i)
    subspace = e_states(E)
    return len(subspace)

def vec2dm(psi):
    """
    Returns density matrix corresponding to state vector psi

    """
    c = psi.column()
    return c * c.conjugate_transpose()

def get_block(A, E):
    """
    Returns E-th block of A in exc basis

    """
    I = J = e_states(E)
    return A[I, J]

def get_full(B, E):
    """
    Returns full matrix A given E-th block of A in exc basis

    Sets all elements of A to zero except elements in block B

    """
    A = matrix(B.base_ring(), states_count)
    I = J = e_states(E)
    A[I, J] = B
    return A

def basis_state(state):
    """
    Returns basis state number ``state``

    """
    psi = vector(CDF, states_count)
    psi[state] = 1
    return psi

def bits(n, min_width = 1):
    """
    Returns binary representation of n as a list of ones and zeroes

    The returned list starts with a major binary digit

    """
    l = Integer(n).bits()
    l.reverse()
    padded_list = [0] * (min_width - len(l)) + l
    return padded_list

def html_var(variable, value):
    """
    Prints variable equals value with formatting

    Works with notebook interface

    """
    html('${} = {}$'.format(variable, latex(value)))

def html_vars(l):
    """
    Prints variables and values from a list of pairs

    Works with notebook interface

    """
    for x in l:
        html_var(*x)

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
