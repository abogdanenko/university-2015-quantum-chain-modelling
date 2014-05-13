load('routines.sage')
load('SymbolicComputationBase.sage')
load('SymbolicComputationPlots.sage')
load('SymbolicComputationInteractive.sage')
load('NumericalParams.sage')
load('NumericalComputationBase.sage')
load('NumericalComputationPlots.sage')
load('NumericalComputationInteractive.sage')
load('NumericalComputationAnimation.sage')
load('NumericalComputation.sage')

qubits_count = 4
states_count = 2^qubits_count
exc_list = range(qubits_count + 1)
states_list = range(states_count)

block_sizes = [binomial(qubits_count, E) for E in exc_list]
left_indices = [sum([block_sizes[k] for k in range(E)])
    for E in exc_list]
right_indices = [left_indices[E] + block_sizes[E]
    for E in exc_list]
# first state with excitation number E
first_indices = [2 ** E - 1 for E in exc_list]
exc_number_rainbow = rainbow(qubits_count + 1)

# alias
SymbolicComputation = SymbolicComputationInteractive
