<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;

class SwapGate extends AbstractTwoQubitGate
{
    public const NAME = 'swap';

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
        return $this->drawTwoQubit($qubits, '×', '×');
    }

    protected function buildMatrix(): NDArray
    {
        return $this->la->array([
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 1, 0, 0],
            [0, 0, 0, 1],
        ], dtype: NDArray::complex64);
    }
}
