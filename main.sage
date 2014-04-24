load('routines.sage')
load('SymbolicComputation.sage')
load('SymbolicComputationInteractive.sage')
load('NumericalParams.sage')
load('NumericalComputation.sage')
load('NumericalComputationInteractive.sage')

qubits_count = 4
states_count = 2^qubits_count
energy_list = range(qubits_count + 1)
states_list = range(states_count)

block_sizes = [binomial(qubits_count, E) for E in energy_list]
left_indices = [sum([block_sizes[k] for k in range(E)])
    for E in energy_list]
right_indices = [left_indices[E] + block_sizes[E]
    for E in energy_list]
energy_rainbow = rainbow(qubits_count + 1)

