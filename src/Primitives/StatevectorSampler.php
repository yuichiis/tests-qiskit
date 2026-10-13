<?php
declare(strict_types=1);

namespace Rindow\Qubit\Primitives;

use Rindow\Qubit\QuantumCircuit;

/**
 * Qiskit互換プリミティブ: qiskit.primitives.StatevectorSampler のPHP版
 *
 * ```php
 * $sampler = new StatevectorSampler();
 * $result = $sampler->run([$qc], shots: 1024)->result();
 * $counts = $result[0]->data->meas->getCounts();
 * ```
 */
class StatevectorSampler
{
    private ?int $defaultSeed;

    public function __construct(?int $seed = null)
    {
        $this->defaultSeed = $seed;
    }

    /**
     * @param array<int,QuantumCircuit> $circuits
     */
    public function run(array $circuits, int $shots = 1024, ?int $seed = null): SamplerJob
    {
        return new SamplerJob($circuits, $shots, $seed ?? $this->defaultSeed);
    }
}
