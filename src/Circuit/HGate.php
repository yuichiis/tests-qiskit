<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;

class HGate extends AbstractSingleQubitGate
{
    public const NAME = 'h';

    public function name(): string
    {
        return self::NAME;
    }

    public function inverse(): self
    {
        return new self(la:$this->la);
    }

    protected function buildMatrix(): NDArray
    {
        $s = 1.0 / sqrt(2.0);
        return $this->la->array(
            [[$s, $s], [$s, -$s]],
            dtype: NDArray::complex64
        );
    }
}
