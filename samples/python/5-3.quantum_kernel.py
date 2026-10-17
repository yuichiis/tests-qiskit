# 量子計算 量子カーネル・量子機械学習（分類の入口）

# 量子回路を特徴写像に使う例題だよ！

# ZZFeatureMap でデータを量子状態に埋め込むよ
# カーネル行列を作って 古典SVMと組み合わせるよ
# qiskit-machine-learning で試せるよ
# 学びポイント：

# データの埋め込み方が性能を決めるよ
# 小さいデータから試そうね

# pip install qiskit-machine-learning scikit-learn が必要だよ
import numpy as np
from qiskit.circuit.library import zz_feature_map
from qiskit_machine_learning.kernels import FidelityQuantumKernel
from qiskit_machine_learning.algorithms import QSVC

# 小さい2クラスのデータ（2次元）を用意するよ
# 対角線で分かれる素直な形にしてあるよ。まずはこれで動きを確かめようね
X = np.array([
    [-1.0, -0.6],
    [-0.8, -1.0],
    [-1.2, -1.1],
    [1.0, 0.6],
    [0.8, 1.0],
    [1.2, 1.1],
])
y = np.array([0, 0, 0, 1, 1, 1])

# 特徴写像：データ1点を 2量子ビットの量子状態に埋め込むよ
# zz_feature_map = 旧ZZFeatureMap だよ
# reps を増やすと表現力が上がるよ（埋め込み方が性能を決めるポイント！）
feature_map = zz_feature_map(feature_dimension=2, reps=2, entanglement="linear")
print("特徴写像の量子ビット数:", feature_map.num_qubits)

# 量子カーネル：2つの量子状態の重なり（内積）をカーネル値にするよ
kernel = FidelityQuantumKernel(feature_map=feature_map)

# 訓練データ同士のカーネル行列を作るよ（サンプル数×サンプル数の対称行列）
K_train = kernel.evaluate(x_vec=X)
print("カーネル行列の形:", K_train.shape)
print(np.round(K_train, 3))

# 量子カーネル + 古典SVM（QSVC）で学習・分類するよ
qsvc = QSVC(quantum_kernel=kernel)
qsvc.fit(X, y)

# 訓練データを分類して正解率を見てみるよ
print(f"訓練データの正解率: {qsvc.score(X, y):.3f}")

# 未知の点も予測してみるよ
X_test = np.array([[-1.0, -0.9], [1.0, 0.9]])
print("予測ラベル:", qsvc.predict(X_test))

# カーネル行列をヒートマップで眺めるよ
import matplotlib.pyplot as plt

plt.imshow(K_train, cmap="viridis")
plt.colorbar(label="Kernel Values")
plt.title("Quantum Kernel Matrix")
plt.show()
