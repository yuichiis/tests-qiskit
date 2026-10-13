# 2. Build your circuit
# Open a Python environment, then run this code to build a Bell state (two entangled qubits).


from qiskit import QuantumCircuit
from qiskit.primitives import StatevectorSampler

print("hello")

qc = QuantumCircuit(2)
qc.h(0)
qc.cx(0, 1)
qc.measure_all()

sampler = StatevectorSampler()
result = sampler.run([qc], shots=1024).result()
print(result[0].data.meas.get_counts())


# {'00': 486, '11': 538}
#
# The expected output is a near-even split between '00' and '11'.



# 3. Visualize your results
# To get a histogram of your results, add the following code to your program.
#
# Uncomment lines 2 and 8 if you are not using Python in a Jupyter notebook
import matplotlib.pyplot as plt
from qiskit.visualization import plot_histogram

counts = result[0].data.meas.get_counts()
plot_histogram(counts)

plt.show()



# 4. See what happens
# Try changing the code to see how it affects the results. For example:
# 
# Add a third qubit by changing to QuantumCircuit(3), and add a second CX gate with qc.cx(1,2). The measurements should then change to 000 and 111, which means all three of these qubits have been entangled.
# 
# See your results shift by adding qc.x(1) to the end of the circuit.

