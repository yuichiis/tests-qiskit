# バーンスタイン・ヴァジラーニのアルゴリズム
# 隠れたビット列 s を 1回で見つける例題だよ！

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

# 隠れた文字列 s = 11 の例
s = "11"
n = len(s)
qc = QuantumCircuit(n + 1, n)

qc.x(n)
qc.h(range(n + 1))

# オラクル：sに応じてCXを貼る
for i, bit in enumerate(reversed(s)):
    if bit == "1":
        qc.cx(i, n)

qc.h(range(n))
qc.measure(range(n), range(n))
print(qc.draw())

sim = AerSimulator()
counts = sim.run(transpile(qc, sim), shots=500).result().get_counts()
print(counts)
# 期待：{'11': 500} みたいに隠れた文字列がそのまま出る
