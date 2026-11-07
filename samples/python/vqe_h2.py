"""
VQE による H2 分子の基底状態エネルギー計算 (Qiskit 2.1 以降)

- ハミルトニアン: H2, STO-3G, 結合長 0.735 Å を Jordan-Wigner + 2量子ビット縮約したもの
  (PySCF なしで動くよう係数を直書きしています)
- ansatz: real_amplitudes (RY + CX)
- 期待値計算: StatevectorEstimator (厳密な状態ベクトル、ショットノイズなし)
- オプティマイザ: SciPy の COBYLA

インストール:  pip install qiskit scipy numpy
"""

import numpy as np
from scipy.optimize import minimize

from qiskit.circuit.library import real_amplitudes
from qiskit.primitives import StatevectorEstimator
from qiskit.quantum_info import SparsePauliOp

# ---------------------------------------------------------------------------
# 1. ハミルトニアン (電子部分のみ)。核間反発エネルギーは別に足す。
# ---------------------------------------------------------------------------
NUCLEAR_REPULSION = 0.7199689944489797  # Hartree (0.735 Å)

hamiltonian = SparsePauliOp.from_list(
    [
        ("II", -1.052373245772859),
        ("IZ", 0.39793742484318045),
        ("ZI", -0.39793742484318045),
        ("ZZ", -0.01128010425623538),
        ("XX", 0.18093119978423156),
    ]
)

# ---------------------------------------------------------------------------
# 2. ansatz (パラメータ付き回路)
# ---------------------------------------------------------------------------
ansatz = real_amplitudes(num_qubits=2, reps=2, entanglement="linear")
print(ansatz.decompose().draw(fold=-1))
print(f"パラメータ数: {ansatz.num_parameters}\n")

# ---------------------------------------------------------------------------
# 3. コスト関数: <ψ(θ)|H|ψ(θ)>
# ---------------------------------------------------------------------------
estimator = StatevectorEstimator()
history = []


def cost(theta: np.ndarray) -> float:
    # Estimator の入力は (回路, 観測量, パラメータ値) の PUB
    pub = (ansatz, hamiltonian, theta)
    result = estimator.run([pub]).result()[0]
    energy = float(result.data.evs)
    history.append(energy)
    return energy


# ---------------------------------------------------------------------------
# 4. 古典オプティマイザで最小化
# ---------------------------------------------------------------------------
rng = np.random.default_rng(seed=42)
theta0 = rng.uniform(-np.pi, np.pi, size=ansatz.num_parameters)

res = minimize(cost, theta0, method="COBYLA", options={"maxiter": 500, "tol": 1e-8})

vqe_electronic = res.fun
vqe_total = vqe_electronic + NUCLEAR_REPULSION

# ---------------------------------------------------------------------------
# 5. 厳密対角化との比較 (2量子ビットなので行列を直接作れる)
# ---------------------------------------------------------------------------
exact_electronic = float(np.linalg.eigvalsh(hamiltonian.to_matrix()).min())
exact_total = exact_electronic + NUCLEAR_REPULSION

print(f"最適化の反復回数 (関数評価回数): {len(history)}")
print(f"VQE   電子エネルギー : {vqe_electronic:.8f} Ha")
print(f"厳密  電子エネルギー : {exact_electronic:.8f} Ha")
print(f"VQE   全エネルギー   : {vqe_total:.8f} Ha")
print(f"厳密  全エネルギー   : {exact_total:.8f} Ha")
print(f"誤差                 : {abs(vqe_total - exact_total):.2e} Ha")