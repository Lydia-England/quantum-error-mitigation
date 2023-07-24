import random 
import numpy as np
import matplotlib as mpl

# from qiskit.quantum_info import DensityMatrix
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise import depolarizing_error
from qiskit import (QuantumCircuit, ClassicalRegister, QuantumRegister, transpile, BasicAer)
from qiskit.circuit.library import (IGate, XGate, YGate, ZGate, HGate, TGate, SGate, CXGate)
from qiskit.visualization import circuit_drawer
from qiskit.result import Result
from qiskit.quantum_info import Statevector

#==== Simulation Parameters ====#
d         = 20     # number of qubits
n         = 6     # circuit depth
epsilon   = 0.01  # error rate
M         = 4000  # total number of runs
gamma_b   = (((1+epsilon/2)/(1-epsilon))**(n*d/2)*(((1+7*epsilon/8)/(1-epsilon))**(n*d/4)))  # simulation overhead
print('gamma_b = ', gamma_b)

#==== Instances of Gates ====#
I        = IGate()   # Identity gate
I.label  = r'I'
X        = IGate()   # X gate
X.label  = r'X'
Y        = IGate()   # Y gate
Y.label  = r'Y'
Z        = IGate()   # Z gate
Z.label  = r'Z'
H        = HGate()   # Hadamard gate
H.label  = r'H'
S        = SGate()   # S gate (Z**0.5)
S.label  = r'S'
T        = TGate()   # T gate (Z**0.25)
T.label  = r'T'
CX       = CXGate()  # CNOT (controlled-X) gate
CX.label = r'CNOT'

#==== List of Gates (to build ideal circuit) ====#
gate_list = [I, H, S, T]
gate_list_names = ['id', 'h', 's', 't']
paulis    = [I, X, Y, Z]



def get_qc_ub(n, d, single_gates):
    qr = QuantumRegister(n, 'qr')  # quantum register 
    qc = QuantumCircuit(qr)        # quantum circuit

    qubit_index = list(range(n))   # list from 0:n
    
    for layer in range(d):
        # odd layer: single qubit gates
        if layer%2 != 0: 
            for qubit in range(n):
                qc.append(random.choice(single_gates), [qubit])
            qc.barrier()
        # even layer: cnot gates
        else:  
            rand_qubit_index = random.shuffle(qubit_index)
            ctrl = qubit_index[:int(n/2)]
            test = qubit_index[int(n/2):]
            qc.cx(ctrl, test)
            qc.barrier()
    return qc

qc = get_qc_ub(n, d, gate_list)
qc.draw('mpl')


# set initial state of the simulator to the ground state using from_int
state = Statevector.from_int(0, 2**n)

# evolve the state by the quantum circuit
state = state.evolve(qc)

# draw using latex
state.draw('latex')


# Probabilities for measuring both qubits
probs = state.probabilities()
print('probs: {}'.format(probs))

# Probabilities for measuring only qubit-0
probs_qubit_0 = state.probabilities([0])
print('Qubit-0 probs: {}'.format(probs_qubit_0))

# Probabilities for measuring only qubit-1
probs_qubit_1 = state.probabilities([1])
print('Qubit-1 probs: {}'.format(probs_qubit_1))

# Probabilities for measuring only qubit-2
probs_qubit_2 = state.probabilities([2])
print('Qubit-2 probs: {}'.format(probs_qubit_2))

# Probabilities for measuring only qubit-3
probs_qubit_3 = state.probabilities([3])
print('Qubit-3 probs: {}'.format(probs_qubit_3))

# Probabilities for measuring only qubit-4
probs_qubit_4 = state.probabilities([4])
print('Qubit-4 probs: {}'.format(probs_qubit_4))

# Probabilities for measuring only qubit-5
probs_qubit_5 = state.probabilities([5])
print('Qubit-5 probs: {}'.format(probs_qubit_5))


#==== Depolarizing Errors ====#
single_err = depolarizing_error(epsilon, 1)
multi_err  = depolarizing_error(epsilon, 2)

#==== Build Noise Model ====#
dk_noise_model = NoiseModel()
dk_noise_model.add_all_qubit_quantum_error(single_err, gate_list_names)
dk_noise_model.add_all_qubit_quantum_error(multi_err, ['CNOT'])
