# 量子フーリエ変換（2量子ビット）QFT: Quantum Fourier Transform

# 位相推定やショア法の部品だよ。まずは2量子ビットで回してみよう！

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit_aer.library import SaveStatevector

qc = QuantumCircuit(2)
qc.h(1)
qc.cp(3.14159 / 2, 0, 1)
qc.h(0)
qc.swap(0, 1)
qc.append(SaveStatevector(2), [0, 1])
print(qc.draw())

sim = AerSimulator(method="statevector")
state = sim.run(transpile(qc, sim)).result().get_statevector()
print(state)
