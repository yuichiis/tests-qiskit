# 量子版コイン投げ
# 2量子ビットのもつれだよ。測ると 00 か 11 しか出ないのが不思議だね

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

qc = QuantumCircuit(1, 1)
qc.h(0)
qc.measure(0, 0)
print(qc.draw())

sim = AerSimulator()
counts = sim.run(transpile(qc, sim), shots=1000).result().get_counts()
print(counts)
# 期待：{'0': 約500, '1': 約500}

