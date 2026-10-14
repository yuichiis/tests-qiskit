# 量子計算 QAOA（Max-Cut・4点リング）

# グラフを2つに分けて 切れる辺を最大化する例題だよ。組み合わせ最適化の入口だよ！

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

import numpy as np

# 4点リング：辺 (0-1, 1-2, 2-3, 3-0)
p = 1
beta, gamma = 0.5, 0.5
qc = QuantumCircuit(4)
qc.h(range(4))

# 問題ユニタリ（ZZ相互作用）
for i, j in [(0, 1), (1, 2), (2, 3), (3, 0)]:
    qc.cx(i, j)
    qc.rz(2 * gamma, j)
    qc.cx(i, j)

# ミキサー
qc.rx(2 * beta, range(4))
qc.measure_all()
print(qc.draw())

sim = AerSimulator()
counts = sim.run(transpile(qc, sim), shots=1000).result().get_counts()
print(counts)

import matplotlib.pyplot as plt
from qiskit.visualization import plot_histogram

plot_histogram(counts)
plt.show()
