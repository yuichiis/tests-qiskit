# グローバー探索（4件から1件探し）MCXGate 版

# 4つの中から当たりを効率よく探す例題だよ。2量子ビットで手計算できるのが嬉しいね
# MCXGate を使うと、複数の制御ビットを持つ X ゲートをまとめて書けるよ！

from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import MCXGate
from qiskit_aer import AerSimulator


# 複数制御の Z ゲート（位相反転）を MCXGate で作る
def mcz(qc, qubits):
    qc.h(qubits[-1])
    qc.append(MCXGate(len(qubits) - 1), qubits)
    qc.h(qubits[-1])


# 当たりを |11> にする例
qc = QuantumCircuit(2, 2)
qc.h([0, 1])

# オラクル：|11> に印をつける
mcz(qc, [0, 1])

# 拡散操作
qc.h([0, 1])
qc.x([0, 1])
mcz(qc, [0, 1])
qc.x([0, 1])
qc.h([0, 1])

qc.measure([0, 1], [0, 1])
print(qc.draw())

sim = AerSimulator()
counts = sim.run(transpile(qc, sim), shots=500).result().get_counts()
print(counts)
# 期待：{'11': ほぼ500} になる
