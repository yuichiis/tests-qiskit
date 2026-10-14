# ドイチュ・ジョサのアルゴリズム
# 関数の性質を 1回の質問で見分ける例題だよ

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

# 平衡関数の例（CXで作る）
qc = QuantumCircuit(2, 1)
qc.h([0, 1])
qc.cx(0, 1)
qc.h(0)
qc.measure(0, 0)
print(qc.draw())

sim = AerSimulator()
counts = sim.run(transpile(qc, sim), shots=500).result().get_counts()
print(counts)
# 期待：平衡関数なら `1` が出る、定数関数なら `0` が出る
