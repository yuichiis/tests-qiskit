<?php
declare(strict_types=1);

namespace Rindow\Qubit\Primitives;

use Rindow\Qubit\QuantumInfo\Statevector;

/**
 * 1回路分(PUB)の実行結果. qiskitの `result[0]` に対応.
 */
class PubResult
{
    public MeasData $meas;
    public Statevector $statevector;

    /**
     * @param array<string,int> $counts
     */
    public function __construct(array $counts, int $numShots, Statevector $statevector)
    {
        $this->meas = new MeasData($counts, $numShots);
        $this->statevector = $statevector;
    }

    /** @return array<string,int> */
    public function getCounts(): array
    {
        return $this->meas->getCounts();
    }

    public function getStatevector(): Statevector
    {
        return $this->statevector;
    }
}
