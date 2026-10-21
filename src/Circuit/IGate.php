<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;

class IGate extends SingleQubitGate
{
    public const NAME = 'id';

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
        return $this->mo->array([[1, 0], [0, 1]], dtype: NDArray::complex128);
    }
}
