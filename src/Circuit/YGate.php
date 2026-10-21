<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;
use function Rindow\Math\Matrix\C;

class YGate extends SingleQubitGate
{
    public const NAME = 'y';

    public function name(): string
    {
        return self::NAME;
    }

    public function inverse(): self
    {
        return new self($this->mo);
    }

    protected function buildMatrix(): NDArray
    {
        return $this->mo->array(
            [[0, C(0.0, -1.0)], [C(0.0, 1.0), 0]],
            dtype: NDArray::complex128
        );
    }
}
