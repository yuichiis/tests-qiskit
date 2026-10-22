<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;
use function Rindow\Math\Matrix\C;

class SdgGate extends AbstractSingleQubitGate
{
    public const NAME = 'sdg';

    public function name(): string
    {
        return self::NAME;
    }

    public function inverse(): SGate
    {
        return new SGate($this->mo);
    }

    protected function buildMatrix(): NDArray
    {
        return $this->mo->array([[1, 0], [0, C(0.0, -1.0)]], dtype: NDArray::complex128);
    }
}
