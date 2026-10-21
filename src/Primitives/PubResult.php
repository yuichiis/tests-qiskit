<?php
declare(strict_types=1);

namespace Rindow\Qubit\Primitives;

use Rindow\Qubit\QuantumInfo\Statevector;

/**
 * 1回路分(PUB)の実行結果. qiskitの `result[0].data.meas` に対応.
 */
class PubResult
{
    public PubData $data;

    /**
     * @param array<string,int> $counts
     */
    public function __construct(array $counts, int $numShots, Statevector $statevector)
    {
        $this->data = new PubData(
            new MeasData($counts, $numShots),
            $statevector
        );
    }
}
