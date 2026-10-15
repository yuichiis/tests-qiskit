# 量子計算 QAOA（Max-Cut・4点リング）

# グラフを2つに分けて 切れる辺を最大化する例題だよ。組み合わせ最適化の入口だよ！
# QAOAは量子回路のパラメータ(beta, gamma)を古典オプティマイザで
# 調整してカット数の期待値を最大化するハイブリッド算法だよ！

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.quantum_info import Statevector
from scipy.optimize import minimize

import numpy as np

# 4点リング：辺 (0-1, 1-2, 2-3, 3-0)
edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
p = 1


def build_qaoa(gamma, beta):
    qc = QuantumCircuit(4)
    qc.h(range(4))

    # 問題ユニタリ（ZZ相互作用）
    for i, j in edges:
        qc.cx(i, j)
        qc.rz(2 * gamma, j)
        qc.cx(i, j)

    # ミキサー
    qc.rx(2 * beta, range(4))
    return qc


def cut_value(bits):
    # qiskitのビット列は左が上位qubitなので右端が q0
    b = [int(bits[-1 - i]) for i in range(4)]
    return sum(1 for i, j in edges if b[i] != b[j])


def neg_expected_cut(params):
    # 状態ベクトルからカット数の期待値を正確に求める (最適化用)
    gamma, beta = params
    probs = Statevector(build_qaoa(gamma, beta)).probabilities_dict()
    return -sum(prob * cut_value(bits) for bits, prob in probs.items())


# 古典オプティマイザでパラメータを最適化 (期待カット数最大化)
res = minimize(neg_expected_cut, x0=[0.5, 0.5], method="COBYLA",
               options={"maxiter": 200})
gamma, beta = res.x
print(f"最適パラメータ: gamma={gamma:.4f}, beta={beta:.4f} (期待カット数 {-res.fun:.4f})")

# 最適パラメータで測定
qc = build_qaoa(gamma, beta)
qc.measure_all()
print(qc.draw())

sim = AerSimulator()
counts = sim.run(transpile(qc, sim), shots=1000).result().get_counts()
print(counts)

import matplotlib.pyplot as plt
from qiskit.visualization import plot_histogram

plot_histogram(counts)
plt.show()

# 答えの表示：カット数（両端が別グループの辺の数）が最大のビット列を探す
best = max(cut_value(bits) for bits in counts)
winners = sorted(bits for bits in counts if cut_value(bits) == best)
print(f"答え: {winners} (カット数 {best} / 全4辺)")
for bits in winners:
    print(f"  {bits}: {counts[bits]}回")
