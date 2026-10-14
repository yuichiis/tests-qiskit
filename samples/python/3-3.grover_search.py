# グローバー探索（4件から1件探し）

# 4つの中から当たりを効率よく探す例題だよ。2量子ビットで手計算できるのが嬉しいね

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

# 当たりを |11> にする例
qc = QuantumCircuit(2, 2)
qc.h([0, 1])

# オラクル：|11> に印をつける
qc.cz(0, 1)

# 拡散操作
qc.h([0, 1])
qc.x([0, 1])
qc.cz(0, 1)
qc.x([0, 1])
qc.h([0, 1])

qc.measure([0, 1], [0, 1])
print(qc.draw())

sim = AerSimulator()
counts = sim.run(transpile(qc, sim), shots=500).result().get_counts()
print(counts)
# 期待：{'11': ほぼ500} になる
