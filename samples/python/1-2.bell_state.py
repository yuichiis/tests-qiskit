# ベル状態
# 2量子ビットのもつれだよ。測ると 00 か 11 しか出ないのが不思議だね

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

qc = QuantumCircuit(2, 2)
qc.h(0)
qc.cx(0, 1)
qc.measure([0, 1], [0, 1])
print(qc.draw())

sim = AerSimulator()
counts = sim.run(transpile(qc, sim), shots=1000).result().get_counts()
print(counts)
# 期待：{'00': 約500, '11': 約500}、`01` と `10` はほぼ0
