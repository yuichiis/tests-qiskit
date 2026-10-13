<?php
declare(strict_types=1);

namespace Rindow\Qubit\Primitives;

/**
 * qiskitの DataBin に相当. 古典レジスタは単一なので meas のみを持つ.
 */
class PubData
{
    public MeasData $meas;

    public function __construct(MeasData $meas)
    {
        $this->meas = $meas;
    }
}
