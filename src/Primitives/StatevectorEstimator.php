<?php
declare(strict_types=1);

namespace Rindow\Qubit\Primitives;

/**
 * qiskit.primitives.StatevectorEstimator のPHP版.
 *
 * ショットノイズのない厳密な状態ベクトルで <ψ|H|ψ> を返す.
 *
 * ```php
 * $estimator = new StatevectorEstimator();
 * $result = $estimator->run([[$ansatz, $hamiltonian, $theta]])->result();
 * $energy = $result[0]->data->evs;
 * ```
 */
class StatevectorEstimator
{
    /**
     * @param array<int,mixed> $pubs PUB の配列. 各PUBは
     *        [QuantumCircuit|RealAmplitudes, SparsePauliOp, array<float>|null]
     */
    public function run(array $pubs): EstimatorJob
    {
        return new EstimatorJob($pubs);
    }
}
