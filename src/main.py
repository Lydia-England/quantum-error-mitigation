#################
###  IMPORTS  ###
#################

# Import My Modules
from twirl import load_pauli_twirling_sets as load_twrl
from twirl import pauli_twirling as twrl
from data  import twirling_groups

# Import General Modules
import datetime as dt
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from math import pi
from itertools import combinations

# Import Qiskit Modules
### importing basic plot tools
from qiskit.visualization       import plot_histogram
from qiskit.tools.visualization import circuit_drawer
### importing Qiskit
from qiskit import (BasicAer, QuantumCircuit, ClassicalRegister, QuantumRegister, transpile)
from qiskit import *
### Operator class
from qiskit.quantum_info import Operator, state_fidelity
### Classes for building up a directed-acyclic graph (DAG) structure
from qiskit.circuit    import Gate
from qiskit.dagcircuit import DAGCircuit
### need gate classes for generating the Pauli twirling sets
from qiskit.circuit.library import (IGate, XGate, YGate, ZGate, CXGate, CZGate, RXGate, RZZGate, SdgGate, ECRGate, iSwapGate) 
### Transpiler stuff needed to make a pass and passmanager
from qiskit.transpiler            import PassManager
from qiskit.transpiler.basepasses import TransformationPass
from qiskit.transpiler.passes     import Optimize1qGatesDecomposition
### A fake system to transpile against
from qiskit.providers.fake_provider import FakeHanoiV2
# import for building ZZ gate
from qiskit.extensions   import  UnitaryGate
backend = BasicAer.get_backend('unitary_simulator')
# import for building noise model
from qiskit.quantum_info import  Kraus, SuperOp
from qiskit_aer          import  AerSimulator
from qiskit_aer.noise    import (NoiseModel, QuantumError, ReadoutError, pauli_error, depolarizing_error, thermal_relaxation_error)


#####################
###  DEFINITIONS  ###
#####################

#--- Define Angles ---#
theta_h = 0
theta_J = -pi/2

#--- Number of Qubits ---#
n = 14

#--- Define Gates ---#
RX_h  = RXGate(theta_h)    # RX(theta h) Gate
RX_h.label  = r'$RX_h$'
SY    = YGate().power(1/2) # Square Root of Y Gate
SY.label    = r'$\sqrt{Y}$'
SYDG  = SY.adjoint()       # SY Gate Adjoint
SYDG.label  = r'$\sqrt{Y}^\dagger$'
SDG   = SdgGate()          # S Gate Adjoint


