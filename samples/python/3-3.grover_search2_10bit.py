# グローバー探索（1024件から1件探し）MCXGate 版

# 10量子ビット（1024通り）の中から当たりを効率よく探す例題だよ
# 当たりが |1111111111> になるようにオラクルを組むよ
# 最適な繰り返し回数は floor(pi/4 * sqrt(1024)) = 25 回！

from math import floor, pi, sqrt

from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import MCXGate
from qiskit_aer import AerSimulator

N = 10
TARGET = '1' * N
ITERATIONS = floor(pi / 4 * sqrt(2 ** N))


# 複数制御の Z ゲート（位相反転）を MCXGate で作る
def mcz(qc, qubits):
    qc.h(qubits[-1])
    qc.append(MCXGate(len(qubits) - 1), qubits)
    qc.h(qubits[-1])


# 当たりを |111...1> にする例
qc = QuantumCircuit(N, N)
qc.h(range(N))

# オラクルと拡散操作を最適な回数だけ繰り返す
for _ in range(ITERATIONS):
    # オラクル：|111...1> に印をつける
    mcz(qc, list(range(N)))

    # 拡散操作
    qc.h(range(N))
    qc.x(range(N))
    mcz(qc, list(range(N)))
    qc.x(range(N))
    qc.h(range(N))

qc.measure(range(N), range(N))
print(qc.draw())
print(f'繰り返し回数: {ITERATIONS}')

sim = AerSimulator()
counts = sim.run(transpile(qc, sim), shots=500).result().get_counts()
print(counts)
# 期待：{'1111111111': ほぼ500} になる
