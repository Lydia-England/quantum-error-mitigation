########################################################################
###  Numerical Simulation of Probabilistic Error Cancellation (PEC)  ###
###  Author:        Lydia England                                    ###
###  Contact:       lydiajoyengland@gmail.com | lydia@ou.edu         ###
###  Creation Date: 07.24.23                                         ###
###  Last Modified: 07.25.23                                         ### 
########################################################################

#########################
###  IMPORT PACKAGES  ###
#########################

import random 
import statistics
import numpy       as     np
import pandas      as     pd
import matplotlib  as     mpl
from   collections import Counter

import qiskit_aer.noise           as      noise
from   qiskit.providers.aer.noise import (NoiseModel, depolarizing_error)
from   qiskit_aer                 import  AerSimulator
from   qiskit.quantum_info        import (Statevector, DensityMatrix)
from   qiskit.circuit.library     import (IGate, XGate, YGate, ZGate, HGate, TGate, SGate, CXGate)
from   qiskit.visualization       import (circuit_drawer, plot_histogram)
from   qiskit.result              import  Result
from   qiskit                     import (QuantumCircuit, ClassicalRegister, QuantumRegister, transpile, BasicAer, Aer, execute)


##############################
###  FUNCTION DEFINITIONS  ###
##############################

#==== Define function to calculate simulation overhead ====#
def calc_sim_overhead(n, d, epsilon):
    """
    Define a function to calculate simulation overhead (referred to as gamma_{beta}).
    Simulation overhead depends on:
        - n        (number of qubits)
        - d        (circuit depth; that is, number of layers of the circuit)
        - epsilon  (error rate)
    The numerators, the denominator, and the exponents are calculated individually.
    The overhead of single qubit gates and for cnots are calculated individually.
    The overall simulation overhead is calculated from these components.
    The output of the function is the simulation overhead, gamma_b.
    """
    denom     = 1 - epsilon                # calc   denominator (same for both terms) 
    num1      = 1 + ( epsilon / 2 )        # calc   numerator of 1-qubit overhead term
    num2      = 1 + ( (7 * epsilon) / 8 )  # calc   numerator of CNOT    overhead term
    exp1      = (n * d) / 2                # calc   exponent  of 1-qubit term (total number of 1-qubit gates)
    exp2      = (n * d) / 4                # calc   exponent  of CNOT    term (total number of CNOT    gates)
    single_gb = (num1 / denom) ** exp1     # calc   simulation overhead for all 1-qubit gates
    cnot_gb   = (num2 / denom) ** exp2     # calc   simulation overhead for all CNOT    gates
    gamma_b   = single_gb * cnot_gb        # calc   total simulation overhead of circuit (gamma_{beta})
    return gamma_b                         # return total simulation overhead of circuit (gamma_{beta})


#==== Define function to create the quantum circuit ====#
def get_qc_ub(n, d, single_gates):
    """
    Define a function to create the quantum circuit.
    This function depends on:
        - n             (number of qubits)
        - d             (circuit depth; that is, number of layers of the circuit)
        - single_gates  (set of single-qubit gates to be randomly applied on alternating circuit layers)
    For every odd-numbered layer, random single-qubit gates from the given list are applied to each qubit.
    For every even-numbered layer, CNOTs are applied to pairs of qubits. 
        - Note: Control and target qubits are chosen randomly.
    After each layer, a barrier is applied for readability. 
    Each layer will thus have eiter n single-qubit gates or n/2 CNOTs.
    The output of the function is the resultant quantum circuit.
    """
    qr = QuantumRegister(n, 'qr')   # create quantum register 
    qc = QuantumCircuit(qr)         # create quantum circuit
    qi = list(range(n))             # list from 0:n
    for layer in range(d):          # loop through all circuit layers (for layer in range (layers))
        if layer%2 != 0:            # odd layer: single qubit gates
            for qubit in range(n):  # add random single qubit gates from list to each qubit
                qc.append(random.choice(single_gates), [qubit])  
            qc.barrier()            # add barrier at end of layer 
        else:                       # even layer: cnot gates
            random.shuffle(qi)      # randomly shuffle list from 0:n
            ctrl = qi[:int(n/2)]    # second half of randomly shuffled qubits is control qubits
            test = qi[int(n/2):]    # first  half of randomly shuffled qubits is target  qubits
            qc.cx(ctrl, test)       # apply CNOT (cx) gate to all ctrl-test pairs
            qc.barrier()            # add barrier at end of layer 
    return qc                       # return the final quantum circuit


#==== Define function to create depolarizing noise model ====#
def get_dk_noise(one_q_gates, two_q_gates, epsilon):
    """
    Define a function to create a depolarizing noise model.
    The depolarizing noise model depends on:
        - one_qubit_gates  (single-qubit gates on which to apply depolarizing noise)
        - two_qubit_gates  (two-qubit    gates on which to apply depolarizing noise)
        - epsilon          (error rate)
    Depolarizing channels are created for the set of single-qubit gates, then for the set of two-qubit gates.
    These channels are added to the noise model (created with NoiseModel() module)
    The output of the function is the whole depolarizing noise model.
    """
    one_err  = depolarizing_error(epsilon, 1)               # create depolarizing channel acting on one-qubit gates
    two_err  = depolarizing_error(epsilon, 2)               # create depolarizing channel acting on two-qubit gates
    dn_mod = NoiseModel()                                   # declare noise model using NoiseModel() module
    dn_mod.add_all_qubit_quantum_error(one_err,one_q_gates) # add one-qubit errors to noise model
    dn_mod.add_all_qubit_quantum_error(two_err,two_q_gates) # add two-qubit errors to noise model
    return   dn_mod                                         # return the final depolarizing noise model


#==== Define function to get projector onto basis states with probabilities above the median value =============#
def get_projector_geq_med(psi):
    """
    Define a function to obtain projector operator projecting onto basis states with probability above the median.
    This projector depends on one argument: psi (statevector).
    The probabilities are extracted from psi; the median value is extracted from psi. 
    A list is created to store 1s and 0s associated with the probabilities for each state.
    If the probability associated with a state is above or equal the median value, a 1 is added to a list.
    If the probability associated with a state is below the median value, a 0 is added to a list. 
    A new statevector is created from the resultant list, and the statevector is converted to an operator.
    The output of this function is that operator, which projects onto the selected subset of basis states.
    """
    probs      = psi.probabilities()                       # get probabilities list from input statevector psi
    median     = statistics.median(probs)                  # identify median value of probabilities list
    state_list = []                                        # initialize list to store 1s and 0s associated w/ basis state probabilities
    for i in range(len(probs)):                            # loop through list of probabilities from psi
        if   probs[i] <  median: state_list.append(0)      # basis state w/ probability less    than     median: add a 0 to the list
        elif probs[i] >= median: state_list.append(1)      # basis state w/ probability greater or equal median: add a 1 to the list
    projector  = Statevector(state_list).to_operator()     # create statevector from state list; convert to projector operator
    return projector                                       # return the resultant projector operator

#==== Define function to randomly generate circuit from given parameters, execute simulations, find delta_0 ====#
def circ_pass(n, d, epsilon, M, gate_list, gate_names):
    #==== Build quantum circuit ==============================#
    """
    Build quantum circuit with function get_qc_ub().
    Quantum circuit has n qubits, d layers.
    Select 1-qubit gates from list: gate_list.
    """
    qc = get_qc_ub(n, d, gate_list)                           # create ideal quantum circuit 
    # qc.draw('mpl')                                            # draw quantum circuit using matplotlib rendering

    #==== Evolve Statevector by quantum circuit ==============#
    psi = Statevector.from_int(0, 2**n)                       # set initial simulator state to ground state using from_int
    psi = psi.evolve(qc)                                      # evolve the state by the quantum circuit
    psi.draw('latex')                                         # draw using latex rendering

    #==== Get projector operator A ===========================#
    """
    Observable A is chosen as a projector onto the subset of 2^(n-1) basis states whose probability
    in the final state of the ideal circuit is above the median value (temme, 2017).
    """
    A = get_projector_geq_med(psi)                             # get projector A from statevector psi

    #==== Calculate ideal expectation values =================#
    # eval_ideal = psi.expectation_value(A)                     # get ideal expectation value from projector A w.r.t. statevector psi

    #===== Add statevector 'screenshot' to end of circuit ====#
    """
    This is necessary so that the simulations below are able to extract the statevector.
    """
    qc.save_statevector()                                     # add save statevector checkpoint to end of qc

    #==== Simulate Ideal Circuit =============================#
    """
    When simulating ideal circuits, changing the method between the exact simulation methods: 
    stabilizer, statevector, density_matrix, and matrix_product_state, should not change the simulation result.
    We choose the statevector simulator here, since it doesn't matter which of the above we choose for the ideal 
    circuit, but we want to obtain the statevector result from the noisy circuit. 
    """
    sim_ideal    = AerSimulator(method='statevector')         # get statevector simulation method backend
    qc_ideal     = transpile(qc, sim_ideal)                   # transpile ideal circuit
    result_ideal = sim_ideal.run(qc_ideal, shots=M).result()  # run ideal simulation 

    #==== Calculate ideal expectation values =================#
    psi_ideal    = result_ideal.get_statevector(0)            # get statevector from noisy circuit
    eval_ideal   = psi_ideal.expectation_value(A)             # get noisy expectation value from projector A w.r.t. noisy psi

    #==== Build depolarizing noise model =====================#
    dk_noise_model = get_dk_noise(gate_names,'cx',epsilon)    # get depolarizing noise model 
    dk_basis_gates = dk_noise_model.basis_gates               # get basis gates from noise model

    #==== Perform a noise simulation =========================#
    sim_noisy      = AerSimulator(noise_model=dk_noise_model,
                                  basis_gates=dk_basis_gates) # set noise model as depolarizing noise model
    qc_noisy       = transpile(qc, sim_noisy)                 # transpile noisy circuit
    result_noisy   = sim_noisy.run(qc_noisy,shots=M).result() # run noisy simulation

    #==== Calculate noisy expectation values =================#
    psi_dk         = result_noisy.get_statevector(0)          # get statevector from noisy circuit
    eval_dk        = psi_dk.expectation_value(A)              # get noisy expectation value from projector A w.r.t. noisy psi

    #==== Get noisy device simulation precision (delta_0) ====#
    delta_zero = abs(eval_dk - eval_ideal)                    # calculate absval of difference between evals (delta_0)

    return delta_zero                                         # return delta_0 


#==== Define a function to create bins =======================#
def create_bins(lower_bound, width, quantity):
    """ create_bins returns an equal-width (distance) partitioning. 
        It returns an ascending list of tuples, representing the intervals.
        A tuple bins[i], i.e. (bins[i][0], bins[i][1])  with i > 0 
        and i < quantity, satisfies the following conditions:
            (1) bins[i][0] + width == bins[i][1]
            (2) bins[i-1][0] + width == bins[i][0] and
                bins[i-1][1] + width == bins[i][1]
    """
    bins    = []
    low_arr = []
    count   = lower_bound 
    for i in range(quantity):
       low_arr.append(count) 
       bins.append((low_arr[i], low_arr[i]+width))
       count = count + width
    return bins


#==== Define function to find index of interval where value ==#
def find_bin(value, bins):
    """ bins is a list of tuples, like [(0,20), (20, 40), (40, 60)],
        binning returns the smallest index i of bins so that
        bin[i][0] <= value < bin[i][1]
    """
    for i in range(0, len(bins)):
        if bins[i][0] <= value < bins[i][1]:
            return i
    return -1

##########################
###  DEFINE VARIABLES  ###
##########################

#==== Simulation Parameters ==============================#
d        = 20    # number of qubits     (20 in paper)
n        = 6     # circuit depth        (6 in paper)
epsilon  = 0.01  # error rate           (0.01 in paper)
M        = 4000  # total number of runs (4000 in paper)
num_circ = 500   # number of randomly generated circuits to test
gamma_b  = calc_sim_overhead(n, d, epsilon) # simulation overhead

#==== Instances of Gates =================================#
I        = IGate()   # Identity    gate
I.label  = r'I'      # Identity    gate label "I"
X        = IGate()   # X (pauli)   gate
X.label  = r'X'      # X (pauli)   gate label "X"
Y        = IGate()   # Y (pauli)   gate
Y.label  = r'Y'      # Y (pauli)   gate label "Y"
Z        = IGate()   # Z (pauli)   gate
Z.label  = r'Z'      # Z (pauli)   gate label "Z"
H        = HGate()   # Hadamard    gate
H.label  = r'H'      # Hadamard    gate label "H"
S        = SGate()   # S (Z**0.5)  gate 
S.label  = r'S'      # S (Z**0.5)  gate label "S"
T        = TGate()   # T (Z**0.25) gate
T.label  = r'T'      # T (Z**0.25) gate label "T"
CX       = CXGate()  # CX (CNOT)   gate
CX.label = r'CNOT'   # CX (CNOT)   gate label "CNOT"

#==== List of Gates (to build ideal circuit) =============#
gate_list    = [I, H, S, T]                               # list of gates which will be randomly added to quantum circuit 
gate_names   = ['id', 'h', 's', 't']                      # list of gate identifiers corresponding to gates in gate_list 
paulis       = [I, X, Y, Z]                               # list of pauli gates 
paulis_names = ['id', 'x', 'y', 'z']                      # list of gate identifiers corresponding to pauli gates

#==== Binning Parameters =================================#
bin_lower_bnd = 0
bin_width     = 0.05
bin_quantity  = 40


###########################
###  EXECUTION SECTION  ###
###########################

#==== Create bins for unmitigated deltas =================#
bins = create_bins(bin_lower_bnd,bin_width,bin_quantity)  # create bins with parameters defined above
bins = pd.IntervalIndex.from_tuples(bins)                 # get pandas IntervalIndex of bins

#==== Run simulations and put deltas into bins ===========#
delta_zero_arr = []                                       # initialize array to hold delta_0 values 
for ii in range(num_circ):                                # perform operation N times; get N delta values
    delta_zero_arr.append(circ_pass(n, d, epsilon, M, gate_list, gate_names))
categorical_delta_zero = pd.cut(delta_zero_arr, bins)
print(categorical_delta_zero)



"""
To do:
    Build estimator for ideal expectation values 
"""


