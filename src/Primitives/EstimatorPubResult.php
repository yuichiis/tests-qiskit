<?php
declare(strict_types=1);

namespace Rindow\Qubit\Primitives;

use Rindow\Qubit\QuantumInfo\Statevector;

/**
 * 1つの PUB の estimator 結果. qiskit の `result[0]->data.evs` に対応.
 */
class EstimatorPubResult
{
    /** evs / stds を持つ qiskit の DataBin 相当 */
    public object $data;
    public Statevector $statevector;

    public function __construct(float $expectationValue, Statevector $statevector, float $std = 0.0)
    {
        $this->data = (object)[
            'evs' => $expectationValue,
            'stds' => $std,
        ];
        $this->statevector = $statevector;
    }

    public function getExpectationValue(): float
    {
        return (float)$this->data->evs;
    }

    public function getStatevector(): Statevector
    {
        return $this->statevector;
    }
}
