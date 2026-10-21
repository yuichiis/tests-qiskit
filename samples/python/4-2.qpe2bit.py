# 2量子ビット位相推定 (QPE: Quantum Phase Estimation)

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.circuit.library import MCPhaseGate, QFTGate

n = 3                # 位相を読み取る精度ビット数 (qubit 0..2)
t0, t1 = n, n + 1    # 対象の2量子ビット (qubit 3, 4)
phi = np.pi / 4      # U は |11> にだけ e^{iφ} を掛ける

qc = QuantumCircuit(n + 2, n)

# 固有状態 |11> を用意 (U|11> = e^{iφ}|11>)
qc.x([t0, t1])

# 補助ビットを重ね合わせ
qc.h(range(n))
qc.barrier()

# 制御ユニタリ: qubit j が U^{2^j} を制御
# 制御2個 (補助ビット j と t0)、ターゲット t1 の MCPhaseGate
for j in range(n):
    qc.append(MCPhaseGate(phi * 2**j, 2), [j, t0, t1])
qc.barrier()

# 逆量子フーリエ変換 (QFTGate を使用)
qc.append(QFTGate(n).inverse(), range(n))
qc.barrier()

qc.measure(range(n), range(n))
print(qc.draw())

sim = AerSimulator()
counts = sim.run(transpile(qc, sim), shots=1024).result().get_counts()
print("測定結果:", counts)
# 期待：φ/2π = 1/8 なので 3ビットで '001' がほぼ 1024
