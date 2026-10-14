# 量子テレポーテーション
# 量子の状態を 別の場所に移す仕組みだよ。通信 + もつれを使うよ！

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

# 3量子ビット：送りたい状態(0) + もつれペア(1,2)
qc = QuantumCircuit(3, 3)

# 送りたい状態を作る（ここでは |1> にしてみる）
qc.x(0)

# もつれペアを作る
qc.h(1)
qc.cx(1, 2)

# 送り側の操作
qc.cx(0, 1)
qc.h(0)
qc.measure(0, 0)
qc.measure(1, 1)

# 受け側の修正（条件付き）
qc.cx(1, 2)
qc.cz(0, 2)
qc.measure(2, 2)

print(qc.draw())
sim = AerSimulator()
counts = sim.run(transpile(qc, sim), shots=500).result().get_counts()
print(counts)
# 期待：受け側(右端)がほぼいつも `1` になる
