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
