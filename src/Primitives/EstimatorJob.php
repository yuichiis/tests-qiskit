<?php
declare(strict_types=1);

namespace Rindow\Qubit\Primitives;

use InvalidArgumentException;
use Rindow\Qubit\Circuit\RealAmplitudes;
use Rindow\Qubit\Quantum\StatevectorSimulator;
use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\QuantumInfo\SparsePauliOp;
use Rindow\Qubit\QuantumInfo\Statevector;

/**
 * sampler.run(...)->result() が返す非同期ジョブの estimator 版.
 * このPHP実装では run() 時点で同期的にシミュレーションを実行する.
 */
class EstimatorJob
{
    /** @var array<int,EstimatorPubResult> */
    private array $pubResults;

    /** @param array<int,mixed> $pubs */
    public function __construct(array $pubs)
    {
        if (count($pubs) === 0) {
            throw new InvalidArgumentException('At least one PUB is required.');
        }
        $this->pubResults = [];
        foreach ($pubs as $pub) {
            $this->pubResults[] = $this->runPub($pub);
        }
    }

    private function runPub(mixed $pub): EstimatorPubResult
    {
        if (!is_array($pub)) {
            throw new InvalidArgumentException('A PUB must be [circuit, observable, parameters?].');
        }
        $circuit = $pub['circuit'] ?? $pub[0] ?? null;
        $observable = $pub['observable'] ?? $pub['observables'] ?? $pub[1] ?? null;
        $parameters = $pub['parameter_values']
            ?? $pub['parameters']
            ?? $pub[2]
            ?? [];

        if (!$observable instanceof SparsePauliOp) {
            throw new InvalidArgumentException('PUB observable must be a SparsePauliOp.');
        }

        if ($circuit instanceof RealAmplitudes) {
            $circuit = $circuit->assignParameters(array_values((array)$parameters));
        }
        if (!$circuit instanceof QuantumCircuit) {
            throw new InvalidArgumentException(
                'PUB circuit must be a QuantumCircuit or a RealAmplitudes ansatz.'
            );
        }
        if ($circuit->numQubits() !== $observable->numQubits()) {
            throw new InvalidArgumentException(
                'Circuit and observable qubit counts do not match.'
            );
        }

        // シミュレータは Precision の既定精度 (既定 float64 + complex128) を使う。
        $sim = new StatevectorSimulator($circuit->numQubits());
        $sim->applyCircuit($circuit);
        $statevector = new Statevector($sim);
        $evs = $observable->expectationValue($statevector);
        return new EstimatorPubResult($evs, $statevector);
    }

    public function result(): EstimatorResult
    {
        return new EstimatorResult($this->pubResults);
    }
}
