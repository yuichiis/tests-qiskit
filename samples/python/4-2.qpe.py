# 位相推定 (QPE: Quantum Phase Estimation)
# T ゲートの位相（45度）を読み取る例題だよ！
# 固有状態 |1> を用意する
# 制御ユニタリを繰り返して 位相を kickback させる
# 逆フーリエ変換 + 測定で位相が分かるよ
# 学びポイント：
# ショア法や化学計算の心臓部だよ
# 量子ビットを増やすと 精度が上がるよ

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
import numpy as np

# 位相推定の回路を作成
n = 3  # 位相を推定する精度 (2^n で割った値が誤差)
qc = QuantumCircuit(n + 1, n)

# 固有状態 |1> を用意 (T|1> = e^{iπ/4}|1>)
qc.x(n)

# 補助ビットを重ね合わせ
qc.h(range(n))
qc.barrier()

# 制御ユニタリ (qubit j が U^{2^j} を制御、U=T)
# T^{2^j} の固有値は e^{iπ2^j/4} なので制御位相でかける
# ここでは T ゲートの位相 φ=π/4 を推定する
for j in range(n):
    qc.cp(np.pi * 2**j / 4, j, n)
qc.barrier()

# 逆量子フーリエ変換 (4-1 の QFT と逆向き・逆符号で、H は下位から)
for j in range(n // 2):
    qc.swap(j, n - 1 - j)
for j in range(n):
    for k in range(j - 1, -1, -1):
        qc.cp(-np.pi / 2**(j - k), k, j)
    qc.h(j)
qc.barrier()

# 測定
qc.measure(range(n), range(n))
print(qc.draw())

# シミュレータで実行
sim = AerSimulator()
transpiled_circuit = transpile(qc, sim)
counts = sim.run(transpiled_circuit, shots=1024).result().get_counts()

print("測定結果:", counts)
