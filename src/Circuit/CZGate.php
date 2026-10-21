<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;

class CZGate extends TwoQubitGate
{
    public const NAME = 'cz';

    public function name(): string
    {
        return self::NAME;
    }

    public function inverse(): self
    {
        return new self($this->mo);
    }

    public function draw(array $qubits): array
    {
        return $this->drawTwoQubit($qubits, '●', '●');
    }

    protected function buildMatrix(): NDArray
    {
        return $this->mo->array([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, -1],
        ], dtype: NDArray::complex128);
    }
}
