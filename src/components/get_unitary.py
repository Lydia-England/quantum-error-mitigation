# Get unitary matrix for a given gate
from qiskit import QuantumCircuit, transpile, BasicAer

def get_unitary(qc, backend=BasicAer.get_backend('unitary_simulator'), nd=3):
    """Recieves a quantum circuit and returns its associated unitary matrix"""
    if not isinstance(qc, QuantumCircuit()):
        raise TypeError('expecting type QuantumCircuit')
    else:
        job = backend.run(transpile(qc, backend))
        return job.result().get_unitary(qc, decimals=nd)

