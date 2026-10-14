# 変分量子固有値法 VQE: Variational Quantum Eigensolver
# 分子の基底エネルギーを求める例題だよ。量子化学の定番だよ

# pip install qiskit-nature pyscf qiskit-algorithms が必要だよ
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver

driver = PySCFDriver(
    atom="H 0 0 0; H 0 0 0.735",
    basis="sto3g",
    charge=0,
    spin=0,
    unit=DistanceUnit.ANGSTROM,
)
problem = driver.run()
print("ハミルトニアンの項数:", len(problem.hamiltonian.second_q_op()))
