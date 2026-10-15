# ショア法（15 = 3 × 5 の因数分解・概念編）
# 本物のショア法は重いので シミュレーションでは考え方だけ掴もうね！
# 15の因数分解は a=7 を選ぶと 7^4 ≡ 1 mod 15 の周期4が見つかるよ
# 周期が分かれば gcd(a^{r/2} ± 1, N) で因数が割れるよ
# シミュレータでは qiskit_algorithms の Factoring を使う方法もあるよ
# 学びポイント：
# 周期探し → 因数分解の流れが大事だよ
# 実機・istsim でも 15や21の小さい数から試そうね

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.circuit.library import UnitaryGate
import numpy as np
from math import gcd
from fractions import Fraction

N = 15
a = 7
n = 3  # カウント用量子ビット数 (周期4は 3ビットでちょうど表せる)
m = 4  # N を表すのに必要なビット数 (2^4 = 16 > 15)
work = list(range(n, n + m))  # 作業用レジスタ [3, 4, 5, 6]

qc = QuantumCircuit(n + m, n)

# 作業レジスタを |1> に (モジュラ指数関数の始点)
qc.x(n)

# カウント側を重ね合わせ
qc.h(range(n))
qc.barrier()

# 制御モジュラ乗算 U^{2^j} (|y> -> |a^{2^j} y mod N>)
# j=0: ×7、j=1: ×7^2 ≡ ×4、j=2: ×7^4 ≡ ×1 (恒等変換なので省略)
# ×4 mod 15 は qubitペアの入れ替えと等価なので Fredkin で書ける
qc.cswap(1, n + 0, n + 2)
qc.cswap(1, n + 1, n + 3)

# ×7 mod 15 は置換ユニタリで作る (|15> は使わないので固定)
perm = np.zeros((16, 16))
for x in range(15):
    perm[(a * x) % N, x] = 1
perm[15, 15] = 1
qc.append(UnitaryGate(perm).control(1), [0] + work)
qc.barrier()

# 逆量子フーリエ変換 (カウント側、4-1 の QFT と逆向き・逆符号で)
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

sim = AerSimulator()
counts = sim.run(transpile(qc, sim), shots=1024).result().get_counts()
print("測定結果:", counts)

# 古典後処理：連分数で周期 r を求めて因数にする
found = set()
for bits in counts:
    r = Fraction(int(bits, 2), 2**n).limit_denominator(N).denominator
    t = r
    while t < N:  # r が真の周期の約数だった場合も倍数で試す
        if t % 2 == 0 and pow(a, t, N) == 1:
            for f in (gcd(pow(a, t // 2) - 1, N), gcd(pow(a, t // 2) + 1, N)):
                if 1 < f < N:
                    found.add(f)
            break
        t += r
print("見つかった因数:", sorted(found))
# 期待：測定は 000/010/100/110 の4通り、因数は [3, 5]

