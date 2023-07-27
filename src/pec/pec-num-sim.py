########################################################################
###  Numerical Simulation of Probabilistic Error Cancellation (PEC)  ###
###  Author:        Lydia England                                    ###
###  Contact:       lydiajoyengland@gmail.com | lydia@ou.edu         ###
###  Creation Date: 07.24.23                                         ###
###  Last Modified: 07.26.23                                         ### 
########################################################################

#########################
###  IMPORT PACKAGES  ###
#########################

import random 
import statistics

import numpy             as      np
import pandas            as      pd
import matplotlib        as      mpl
import matplotlib.pyplot as      plt
from   collections       import  Counter

import qiskit_aer.noise           as      noise
from   qiskit.providers.aer.noise import (NoiseModel, depolarizing_error)
from   qiskit_aer                 import  AerSimulator
from   qiskit.quantum_info        import (Statevector, DensityMatrix)
from   qiskit.circuit.library     import (IGate, XGate, YGate, ZGate, HGate, TGate, SGate, CXGate)
from   qiskit.visualization       import (circuit_drawer, plot_histogram)
from   qiskit.result              import  Result
from   qiskit                     import (QuantumCircuit, ClassicalRegister, QuantumRegister, transpile, 
                                          BasicAer, Aer, execute)

from   functions_pec              import (calc_sim_overhead, get_dk_noise, ideal_sim, noise_sim, calc_delta_zero, 
                                          create_bins, find_bin, get_qc_ub, get_projector_geq_median)

##########################
###  DEFINE VARIABLES  ###
##########################

#==== Simulation Parameters ==============================#
d             = 20                                        # number of qubits     (20 in paper)
n             = 6                                         # circuit depth        (6 in paper)
epsilon       = 0.01                                      # error rate           (0.01 in paper)
M             = 4000                                      # total number of runs (4000 in paper)
num_circ      = 50                                       # number of randomly generated circuits to test
gamma_b       = calc_sim_overhead(n, d, epsilon)          # simulation overhead

#==== Instances of Gates =================================#
I             = IGate()                                   # Identity    gate
I.label       = r'I'                                      # Identity    gate label "I"
X             = IGate()                                   # X (pauli)   gate
X.label       = r'X'                                      # X (pauli)   gate label "X"
Y             = IGate()                                   # Y (pauli)   gate
Y.label       = r'Y'                                      # Y (pauli)   gate label "Y"
Z             = IGate()                                   # Z (pauli)   gate
Z.label       = r'Z'                                      # Z (pauli)   gate label "Z"
H             = HGate()                                   # Hadamard    gate
H.label       = r'H'                                      # Hadamard    gate label "H"
S             = SGate()                                   # S (Z**0.5)  gate 
S.label       = r'S'                                      # S (Z**0.5)  gate label "S"
T             = TGate()                                   # T (Z**0.25) gate
T.label       = r'T'                                      # T (Z**0.25) gate label "T"
CX            = CXGate()                                  # CX (CNOT)   gate
CX.label      = r'CNOT'                                   # CX (CNOT)   gate label "CNOT"

#==== List of Gates (to build ideal circuit) =============#
gate_list     = [I, H, S, T]                              # list of gates which will be randomly added to quantum circuit 
gate_names    = ['id', 'h', 's', 't']                     # list of gate identifiers corresponding to gates in gate_list 
paulis        = [I, X, Y, Z]                              # list of pauli gates 
paulis_names  = ['id', 'x', 'y', 'z']                     # list of gate identifiers corresponding to pauli gates

#==== Binning Parameters =================================#
bin_lower_bnd = 0                                         # lower bound of bins
bin_width     = 0.05                                      # width of each bin
bin_quantity  = 40                                        # quantity of bins (total number of bins)


###########################
###  EXECUTION SECTION  ###
###########################


#==== Create bins for unmitigated deltas =================#
bins = create_bins(bin_lower_bnd,bin_width,bin_quantity)  # create bins with parameters defined above
bins = pd.IntervalIndex.from_tuples(bins)                 # get pandas IntervalIndex of bins

#==== Run simulations and put deltas into bins ===========#
delta_zero_arr = []                                       # initialize array to hold delta_0 values 
eval_init_arr  = []                                       # initialize array to hold first calculated ideal ev 
eval_ideal_arr = []                                       # initialize array to hold ideal ev values
eval_noise_arr = []                                       # initialize array to hold noisy ev values

for ii in range(num_circ):                                # perform operation N times; get N delta values
    #==== Build quantum circuit ==========================#
    """
    Build quantum circuit with n qubits, d layers; use function get_qc_ub().
    Select 1-qubit gates from list: gate_list.
    """
    qc  = get_qc_ub(n, d, gate_list)                      # create ideal quantum circuit 
    # qc.draw('mpl')                                      # draw quantum circuit using matplotlib rendering

    #==== Evolve Statevector by quantum circuit ==========#
    psi = Statevector.from_int(0, 2**n)                   # set initial simulator state to ground state using from_int
    psi = psi.evolve(qc)                                  # evolve the state by the quantum circuit
    psi.draw('latex')                                     # draw using latex rendering

    #==== Get projector operator A =======================#
    """
    Observable A is chosen as a projector onto the subset of 2^(n-1) basis states whose probability
    in the final state of the ideal circuit is above the median value (temme, 2017).
    """
    A   = get_projector_geq_median(psi)                   # get projector A from statevector psi

    #==== Calculate ideal expectation values =============#
    eval_ideal_temp = psi.expectation_value(A)               # get ideal expectation value from projector A w.r.t. statevector psi
    eval_init_arr.append(eval_ideal_temp)

    #===== Add statevector 'screenshot' end of circuit ===#
    """
    This is necessary so that the simulations below are able to extract the statevector.
    """
    qc.save_statevector()                                 # add save statevector checkpoint to end of qc
    #==== Build depolarizing noise model =================#
    dk_model = get_dk_noise(gate_names,'cx',epsilon)      # get depolarizing noise model 
    dk_gates = dk_model.basis_gates                       # get basis gates from depolarizing noise model
    ev_ideal = ideal_sim(qc,A,M)                          # simulate ideal circ;                 M shots; get ev of observ. A
    ev_noise = noise_sim(qc,A,M, dk_model, dk_gates)      # simulate circ w/ depolarizing noise; M shots; get ev of observ. A
    dz       = calc_delta_zero(ev_ideal, ev_noise)        # find  delta_0 for random circuit given defined parameters
    print("delta_0 ", ii, " is: ", dz)                    # print delta_0 value 
    delta_zero_arr.append(dz)                             # add   delta_0 to list of delta_0 values
    eval_ideal_arr.append(ev_ideal)                       # add ideal expextation valye of A to list of ideal ev values
    eval_noise_arr.append(ev_noise)                       # add noisy expectation value of A to list of noisy ev values

categorical_dz = pd.cut(delta_zero_arr, bins)             # cut data into categorical object based on bins
print("categorical dz is: ", categorical_dz)              # print categorical object of deltas and bins
dz_val_counts  = pd.value_counts(categorical_dz)          # get counts within each bin from the categorical object
print("dz val counts is: ",  dz_val_counts)               # print counts of delta_0 in each bin

print("Initial ideal expectation values are:", eval_init_arr )
print("Ideal expectation values are:",         eval_ideal_arr)
print("Noisy expectation values are:",         eval_noise_arr)


eval_init_bins = create_bins(0, int(max(eval_init_arr)), int(max(eval_init_arr))*2)


fig, ax = plt.subplots()




"""
To do:
    Build estimator for ideal expectation values 
"""


