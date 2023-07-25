########################################################################
###  Numerical Simulation of Probabilistic Error Cancellation (PEC)  ###
########################################################################

#########################
###  IMPORT PACKAGES  ###
#########################

import random 
import numpy      as np
import matplotlib as mpl

# from qiskit.quantum_info import DensityMatrix
from qiskit.providers.aer.noise import  NoiseModel
from qiskit.providers.aer.noise import  depolarizing_error
from qiskit.quantum_info        import  Statevector
from qiskit.circuit.library     import (IGate, XGate, YGate, ZGate, HGate, TGate, SGate, CXGate)
from qiskit.visualization       import  circuit_drawer
from qiskit.result              import  Result
from qiskit                     import (QuantumCircuit, ClassicalRegister, QuantumRegister, transpile, BasicAer)



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
    qr = QuantumRegister(n, 'qr')        # create quantum register 
    qc = QuantumCircuit(qr)              # create quantum circuit
    qi = list(range(n))                  # list from 0:n
    for layer in range(d):               # loop through all circuit layers (for layer in range (layers))
        if layer%2 != 0:                 # odd layer: single qubit gates
            for qubit in range(n):       # add random single qubit gates from list to each qubit
                qc.append(random.choice(single_gates), [qubit])  
            qc.barrier()                 # add barrier at end of layer 
        else:                            # even layer: cnot gates
            rand_qi = random.shuffle(qi) # randomly shuffle list from 0:n
            ctrl = rand_qi[:int(n/2) ]   # second half of randomly shuffled qubits is control qubits
            test = rand_qi[ int(n/2):]   # first  half of randomly shuffled qubits is target  qubits
            qc.cx(ctrl, test)            # apply CNOT (cx) gate to all ctrl-test pairs
            qc.barrier()                 # add barrier at end of layer 
    return qc                            # return the final quantum circuit


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
    one_err  = depolarizing_error(epsilon, 1)                  # create depolarizing channel acting on one-qubit gates
    two_err  = depolarizing_error(epsilon, 2)                  # create depolarizing channel acting on two-qubit gates
    dn_model = NoiseModel()                                    # declare noise model using NoiseModel() module
    dn_model.add_all_qubit_quantum_error(one_err, one_q_gates) # add one-qubit errors to noise model
    dn_model.add_all_qubit_quantum_error(two_err, two_q_gates) # add two-qubit errors to noise model
    return   dn_model                                          # return the final depolarizing noise model



##########################
###  DEFINE VARIABLES  ###
##########################

#==== Simulation Parameters ====#
d       = 20    # number of qubits
n       = 6     # circuit depth
epsilon = 0.01  # error rate
M       = 4000  # total number of runs
gamma_b = calc_sim_overhead(n, d, epsilon) # simulation overhead

print('gamma_b = ', gamma_b)   # print \gamma_{\beta}

#==== Instances of Gates ====#
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

#==== List of Gates (to build ideal circuit) ====#
gate_list       = [I, H, S, T]           # list of gates which will be randomly added to quantum circuit 
gate_list_names = ['id', 'h', 's', 't']  # list of gate identifiers corresponding to gates in gate_list 
paulis          = [I, X, Y, Z]           # list of pauli gates 
paulis_names    = ['id', 'x', 'y', 'z']  # list of gate identifiers corresponding to pauli gates


###############################
###  BUILD QUANTUM CIRCUIT  ###
###############################
"""
Build quantum circuit with function get_qc_ub().
Quantum circuit has n qubits, d layers.
Select 1-qubit gates from list: gate_list.
"""
qc = get_qc_ub(n, d, gate_list)  # create ideal quantum circuit 
qc.draw('mpl')                   # draw quantum circuit using matplotlib rendering


###############################################
###  EVOLVE STATEVECTOR BY QUANTUM CIRCUIT  ###
###############################################
psi = Statevector.from_int(0, 2**n)  # set initial simulator state to ground state using from_int
psi = psi.evolve(qc)                 # evolve the state by the quantum circuit
psi.draw('latex')                    # draw using latex


# Probabilities for measuring both qubits
probs = psi.probabilities()
print('probs: {}'.format(probs))

# Probabilities for measuring only qubit-0
probs_qubit_0 = psi.probabilities([0])
print('Qubit-0 probs: {}'.format(probs_qubit_0))

# Probabilities for measuring only qubit-1
probs_qubit_1 = psi.probabilities([1])
print('Qubit-1 probs: {}'.format(probs_qubit_1))

# Probabilities for measuring only qubit-2
probs_qubit_2 = psi.probabilities([2])
print('Qubit-2 probs: {}'.format(probs_qubit_2))

# Probabilities for measuring only qubit-3
probs_qubit_3 = psi.probabilities([3])
print('Qubit-3 probs: {}'.format(probs_qubit_3))

# Probabilities for measuring only qubit-4
probs_qubit_4 = psi.probabilities([4])
print('Qubit-4 probs: {}'.format(probs_qubit_4))

# Probabilities for measuring only qubit-5
probs_qubit_5 = psi.probabilities([5])
print('Qubit-5 probs: {}'.format(probs_qubit_5))



########################################
###  BUILD DEPOLARIZING NOISE MODEL  ###
########################################

dk_noise_model = get_dk_noise(gate_list_names, 'CNOT', epsilon)


