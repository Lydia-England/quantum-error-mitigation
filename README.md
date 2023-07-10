# quantum-error-mitigation
Quantum circuit error mitigation analysis modelled after the IBM eagle chip and the recent Nature article by Kim, et. al. (https://doi.org/10.1038/s41586-023-06096-3)


# Structure of the processor, circuit, etc. 

Let's take a deep dive into the (Kim, et. al.) paper so we start on the same page.
I'm assuming you have *some* background in quantum computation?
If not, check out [Quantum Computation and Quantum Information by Nielsen and Chuang](https://www.cambridge.org/highereducation/books/quantum-computation-and-quantum-information/01E10196D0A682A6AEFFEA52D53BE9AE#overview). 

I'll start by listing some key terms, then we can talk about the quantum gates involved in the circuit. 
This information is largely useless until we talk about the circuit itself, but I want to formalize the "language" we're using right away. 

Next, we'll look at the structure of the quantum processor itself (with some unecessary detail, perhaps, but I hope it's interesting). 
A lot of this background is necessary to understand how we'll mathematically describe the processor.
Eventually, it may be useful to analyze different processor topologies, but for now we'll just look at the 127-qubit IBM chip.

Once we've established some knowledge of the processor, we'll look at the benchmark circuit used: the Trotterized time evolution of a 2D transverse-field Ising model, sharing the topology of the qubit processor (Kim, et. al). 
That is, we'll look at the 2D transverse-field quantum Ising model, and then we'll look at Trotterization.
At this point, the gates I will have told you about at the beginning will all be put into context (you're welcome).

We'll then need to pull back and talk about what quantum errors even *are*, and how do we mathematically describe them or model them?
We're trying to mitigate *noise* in our circuit, so we'll need a model for that. 
The model we'll use is the "sparse Pauli-Linblad noise model."
(Fair warning: it's a bit mathematically complicated - or at least hard to explain).

*Finally*, we'll talk about what we all came here for: quantum error mitigation. 
We'll be looking at "zero-noise extrapolation (ZNE)" which we can use to amplify the noise to different gain levels and then estimate the zero-noise limit using extrapolation.


## Quantum gates 
Quantum gates are the basic circuit elements operating on qubits - analagous to classical logic gates in digital circuits.
Quantum gates are unitary operators; that is, they are linear operators on a Hilbert space that preserve the inner product.
We describe these unitary operators as unitary matrices in some basis, usually the *computational basis*.

There are several essential quantum gates which describe the time-evolution of the 127-qubit superconducting quantum processor in (Kim, et. al).
I'll tell you what these *are* now - these are the Pauli gates (I, X, Y, Z), the S gate, the CNOT gate, and the RZZ gate - but don't worry if you don't exactly *how* these are implemented yet. 

### Single-qubit gates
The primary single-qubit gates to know are the identity gate (I), and the Pauli gates (X,Y,Z). 
These are used in the Ising model and in Pauli Twirling (and in just about everything else if you're talking quantum computation).
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

Additionally, we will use the "S-dagger" gate (adjoint of the S gate):
```math
S^{\dagger} =
\begin{pmatrix}
  1 & 0 \\
  0 & -i
\end{pmatrix}
```

We will use the square root of the Y gate, and the adjoint of the square root of the Y gate(these are used to build the RZZ gate):
```math
\sqrt{Y} = \frac{1}{2}
\begin{pmatrix}
  1+i & -1-i \\
  1+i & 1+i
\end{pmatrix}
```
```math
\sqrt{Y}^{\dagger} =
\begin{pmatrix}
  1-i & 1-i \\
  -1+i & 1-i
\end{pmatrix}
```

Finally, we will use the RX gate. 
This gate represents a rotation around the X-axis (in the Bloch sphere representation).
This gate is part of the Trotterization (think step-wise time-evolution) of the Ising model.
```math
R_X(\theta) = 
\begin{pmatrix}
    \cos{\left(\frac{\theta}{2}\right)} & -i \sin{\left(\frac{\theta}{2}\right)} \\
    -i \sin{\left(\frac{\theta}{2}\right)} & \cos{\left(\frac{\theta}{2}\right)}
\end{pmatrix}
```

### Multi-qubit gates

The first essential multi-qubit gate is the CNOT (controlled NOT) gate. 
The action of the control gate is this: if the control qubit is set to 0, the target qubit remains unchanged. 
If the control qubit is set to q, the target qubit is flipped. 
The CNOT gate is represented as:
```math
CNOT = 
\begin{pmatrix}
  1 & 0 & 0 & 0 \\
  0 & 1 & 0 & 0 \\
  0 & 0 & 0 & 1 \\
  0 & 0 & 1 & 0 \\
\end{pmatrix}
```

Finally, we use the RZZ gate. 
This gate is a two-qubit rotation gate.
Its general form is the following:
```math
R_{ZZ}\left(\theta\right) = e^\left(-i\frac{\theta}{2} Z\otimes Z\right) = 
\begin{pmatrix}
    e^{-i\frac{\theta}{2}} & 0 & 0 & 0 \\
    0 & e^{i\frac{\theta}{2}} & 0 & 0 \\
    0 & 0 & e^{i\frac{\theta}{2}} & 0 \\
    0 & 0 & 0 & e^{-i\frac{\theta}{2}} \\
\end{pmatrix}
```

In (Kim, et. al.), only the case where `$\theta_J = -\pi/2$` is considered. 
Thus, the ZZ rotation can be constructed entirely from single qubit gates and one CNOT gate. 


## Topology of the qubit processor


