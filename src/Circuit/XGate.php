<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;

class XGate extends AbstractSingleQubitGate
{
    public const NAME = 'x';

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
        return $this->la->array([[0, 1], [1, 0]], dtype: NDArray::complex64);
    }
}
