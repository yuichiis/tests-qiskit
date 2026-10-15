# 例題2-2：超密度符号（2ビットを1量子ビットで送る）
# もつれがあると 1つの量子ビットで 2ビット分の情報を送れるよ！
#
# 送りたい2ビットで操作を変えるよ：00→何もしない、01→Z、10→X、11→XZ
# 受け側で CX + H + 測定 すると元に戻るよ
# 学びポイント：
#
# もつれは通信の資源になるよ
# テレポートとちょうど逆の関係だよ


from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

# Alice と Bob がもつれを共有した状態からスタート
qc = QuantumCircuit(2, 2)
qc.h(0)
qc.cx(0, 1)
qc.barrier()

# Alice が送りたいビット (例: 10)
message = "10"

# 00: do nothing
if message == "00":
    pass
# 01: Z
elif message == "01":
    qc.z(0)
# 10: X
elif message == "10":
    qc.x(0)
# 11: XZ
elif message == "11":
    qc.x(0)
    qc.z(0)

qc.barrier()

# Bob 側: CX + H + 測定
qc.cx(0, 1)
qc.h(0)
qc.measure([0, 1], [0, 1])

print(qc.draw())

sim = AerSimulator()
counts = sim.run(transpile(qc, sim), shots=1024).result().get_counts()
print(counts)

import matplotlib.pyplot as plt
from qiskit.visualization import plot_histogram

plot_histogram(counts)
plt.show()
# 期待：送ったメッセージと同じビット列が100%出る（例：10なら {'10': 1024}）
