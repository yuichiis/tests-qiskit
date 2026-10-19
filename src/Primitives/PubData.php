<?php
declare(strict_types=1);

namespace Rindow\Qubit\Primitives;

/**
 * qiskitの DataBin に相当. 古典レジスタは単一なので meas を持つ.
 * 状態ベクトルは statevector から取得できる.
 */
class PubData
{
    public MeasData $meas;
    public StatevectorData $statevector;

    public function __construct(MeasData $meas, StatevectorData $statevector)
    {
        $this->meas = $meas;
        $this->statevector = $statevector;
    }
}
