<?php
declare(strict_types=1);

namespace Rindow\Qubit\Primitives;

use Rindow\Qubit\QuantumInfo\Statevector;

/**
 * qiskitの DataBin に相当. 古典レジスタは単一なので meas を持つ.
 * 状態ベクトルは statevector から取得できる.
 */
class PubData
{
    public MeasData $meas;
    public Statevector $statevector;

    public function __construct(MeasData $meas, Statevector $statevector)
    {
        $this->meas = $meas;
        $this->statevector = $statevector;
    }
}
