<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;
use function Rindow\Math\Matrix\C;

class TGate extends AbstractSingleQubitGate
{
    public const NAME = 't';

    public function name(): string
    {
        return self::NAME;
    }

    public function inverse(): TdgGate
    {
        return new TdgGate($this->mo);
    }

    protected function buildMatrix(): NDArray
    {
        $s = 1.0 / sqrt(2.0);
        return $this->mo->array([[1, 0], [0, C($s, $s)]], dtype: NDArray::complex128);
    }
}
