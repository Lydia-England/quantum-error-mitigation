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


#==== Define function to calculate simulation overhead =====#
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
    denom     = 1 - epsilon                                 # calc   denominator (same for both terms) 
    num1      = 1 + ( epsilon / 2 )                         # calc   numerator of 1-qubit overhead term
    num2      = 1 + ( (7 * epsilon) / 8 )                   # calc   numerator of CNOT    overhead term
    exp1      = (n * d) / 2                                 # calc   exponent  of 1-qubit term (total number of 1-qubit gates)
    exp2      = (n * d) / 4                                 # calc   exponent  of CNOT    term (total number of CNOT    gates)
    single_gb = (num1 / denom) ** exp1                      # calc   simulation overhead for all 1-qubit gates
    cnot_gb   = (num2 / denom) ** exp2                      # calc   simulation overhead for all CNOT    gates
    gamma_b   = single_gb * cnot_gb                         # calc   total simulation overhead of circuit (gamma_{beta})
    return gamma_b                                          # return total simulation overhead of circuit (gamma_{beta})


#==== Define function to create the quantum circuit ========#
def get_qc_ub(n, d, single_gates, drawstyle='none', filename='qc'):
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
    qr = QuantumRegister(n, 'qr')                           # create quantum register 
    qc = QuantumCircuit(qr)                                 # create quantum circuit
    qi = list(range(n))                                     # list from 0:n
    for layer in range(d):                                  # loop through all circuit layers (for layer in (layers))
        if layer%2 != 0:                                    # odd layer: single qubit gates
            for qubit in range(n):                          # add random single qubit gates from list to each qubit
                qc.append(random.choice(single_gates), [qubit])  
            qc.barrier()                                    # add barrier at end of layer 
        else:                                               # even layer: cnot gates
            random.shuffle(qi)                              # randomly shuffle list from 0:n
            ctrl = qi[:int(n/2)]                            # second half of randomly shuffled qubits is ctrl   qubits
            test = qi[int(n/2):]                            # first  half of randomly shuffled qubits is target qubits
            qc.cx(ctrl, test)                               # apply CNOT (cx) gate to all ctrl-test pairs
            qc.barrier()                                    # add barrier at end of layer 
    #==== Save circuit diagram to file if applicable =======#
    if   drawstyle=='latex_source': qc.draw('latex_source', filename='./'+filename+'.tex') 
    elif drawstyle=='mpl':          qc.draw('mpl',          filename='./'+filename+'.png') 
    elif drawstyle=='txt':          qc.draw('text',         filename='./'+filename+'.txt') 
    elif drawstyle=='latex':        qc.draw('latex',        filename='./'+filename+'.pdf') 
    return qc                                               # return the final quantum circuit



#==== Define function to create depolarizing noise model ===#
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
def get_projector_geq_median(psi):
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
    probs      = psi.probabilities()                          # get probabilities list from input statevector psi
    median     = statistics.median(probs)                     # identify median value of probabilities list
    state_list = []                                           # init. list to store 1s and 0s associated w/ basis state probs
    for i in range(len(probs)):                               # loop through list of probabilities from psi
        if   probs[i] <  median: state_list.append(0)         # basis state w/ prob less    than     median: add a 0 to the list
        elif probs[i] >= median: state_list.append(1)         # basis state w/ prob greater or equal median: add a 1 to the list
    projector  = Statevector(state_list).to_operator()        # create statevector from state list; convert to proj operator
    return projector                                          # return the resultant projector operator



#==== Define a function to run ideal sim, get exp. val of A ==#
def ideal_sim(qc, A, M):
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
    return eval_ideal


#==== Define a function to run noisy sim, get exp. val of A ==#
def noise_sim(qc, A, M, noise_model, basis_gates):
    #==== Perform a noise simulation =========================#
    sim_noise    = AerSimulator(noise_model=noise_model,
                                basis_gates=basis_gates)      # set noise model as depolarizing noise model
    qc_noise     = transpile(qc, sim_noise)                   # transpile noisy circuit
    result_noise = sim_noise.run(qc_noise,shots=M).result()   # run noisy simulation
    #==== Calculate noisy expectation values =================#
    psi_noise  = result_noise.get_statevector(0)              # get statevector from noisy circuit
    eval_noise = psi_noise.expectation_value(A)               # get noisy expectation value from projector A w.r.t. noisy psi
    return eval_noise


#==== Define a function to calculate delta_0 =================#
def calc_delta_zero(eval_ideal, eval_noise):
    #==== Get noisy device simulation precision (delta_0) ====#
    delta_zero = abs(eval_noise - eval_ideal)                 # calculate absval of difference between evals (delta_0)
    return delta_zero                                         # return delta_0 


