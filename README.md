# quantum-error-mitigation
Quantum circuit error mitigation analysis modelled after the IBM eagle chip and the recent Nature article by Kim, et. al. (https://doi.org/10.1038/s41586-023-06096-3)

## Quantum gates 
Quantum gates are the basic circuit elements operating on qubits - analagous to classical logic gates in digital circuits.
Quantum gates are unitary operators; that is, they are linear operators on a Hilbert space that preserve the inner product.
We describe these unitary operators as unitary matrices in some basis, usually the *computational basis*.

There are several essential quantum gates which describe the time-evolution of the 127-qubit superconducting quantum processor in (Kim, et. al.).
I'll tell you what these *are* now - these are the Pauli gates (I, X, Y, Z), the S gate, the CNOT gate, and the RZZ gate - but don't worry if you don't exactly *how* these are implemented yet. 

### Single-qubit gates
The single-qubit gates to know are the identity gate (I), and the Pauli gates (X,Y,Z). 
These are as follows (in the computational basis):

```math
I = 
\begin{pmatrix}
  1 & 0 \\
  0 & 1
\end{pmatrix}
```
```math
X = \sigma_x = 
\begin{pmatrix}
  0 & 1 \\
  1 & 0
\end{pmatrix}
```
```math
Y = \sigma_y = 
\begin{pmatrix}
  0 & -i \\
  i & 0
\end{pmatrix}
```
```math
Z = \sigma_z = 
\begin{pmatrix}
  1 & 0 \\
  0 & -1
\end{pmatrix}
```
These matrixes are all involutory. 
The Pauli matrixes anti-commute. 
The matrix exponential of a Pauli matrix is a rotation operator.


### Multi-qubit gates

