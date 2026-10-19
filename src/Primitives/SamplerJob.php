<?php
declare(strict_types=1);

namespace Rindow\Qubit\Primitives;

use InvalidArgumentException;
use Rindow\Qubit\Quantum\StatevectorSimulator;
use Rindow\Qubit\QuantumCircuit;

/**
 * sampler.run(...)->result() が返す非同期ジョブの簡易版.
 * このPHP実装では run() 時点で同期的にシミュレーションを実行する.
 */
class SamplerJob
{
    /** @var array<int,PubResult> */
    private array $pubResults;

    /** @param array<int,QuantumCircuit> $circuits */
    public function __construct(array $circuits, int $shots = 1024, ?int $seed = null)
    {
        if ($shots < 1) {
            throw new InvalidArgumentException('shots must be >= 1.');
        }
        if (count($circuits) === 0) {
            throw new InvalidArgumentException('At least one circuit is required.');
        }
        $this->pubResults = [];
        foreach ($circuits as $i => $circuit) {
            if (!$circuit instanceof QuantumCircuit) {
                throw new InvalidArgumentException('run() expects array of QuantumCircuit.');
            }
            // 回路ごとに異なるseed系列にする (qiskitと同様に再現性を保つ)
            $circuitSeed = $seed === null ? null : $seed + $i;
            $sim = new StatevectorSimulator($circuit->numQubits());
            $sim->applyCircuit($circuit);
            $statevector = $sim->getStatevector();
            $map = $circuit->effectiveMeasureMap();
            $numClbits = $circuit->effectiveNumClbits();
            $counts = $sim->sampleCounts($map, $numClbits, $shots, $circuitSeed);
            $this->pubResults[] = new PubResult($counts, $shots, $statevector);
        }
    }

    public function result(): SamplerResult
    {
        return new SamplerResult($this->pubResults);
    }
}
