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
