chain_len = 2
qubits_count = 2 * chain_len + 1
states_count = 2^qubits_count
exc_list = range(qubits_count + 1)
states_list = range(states_count)

load('routines.sage')
load('Subspace.sage')
load('SymbolicComputationBase.sage')
load('SymbolicComputationInteractive.sage')
load('NumericalParams.sage')
load('MEIntegrator.sage')
load('NumericalComputationBase.sage')
load('NumericalComputationPlots.sage')
load('NumericalComputationInteractive.sage')
load('NumericalComputationAnimation.sage')
load('NumericalComputation.sage')
load('Conductivity.sage')

# alias
SymbolicComputation = SymbolicComputationInteractive
