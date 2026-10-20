# 量子フーリエ変換（2量子ビット）QFT: Quantum Fourier Transform

# 位相推定やショア法の部品だよ。まずは2量子ビットで回してみよう！
# QFTGate を使うと、手書きの H や制御位相ゲートをまとめて書けるよ！

from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import QFTGate
from qiskit_aer import AerSimulator
from qiskit_aer.library import SaveStatevector

qc = QuantumCircuit(2)
qc.append(QFTGate(2), [0, 1])
qc.append(SaveStatevector(2), [0, 1])
print(qc.draw())

sim = AerSimulator(method="statevector")
state = sim.run(transpile(qc, sim)).result().get_statevector()
print(state)
