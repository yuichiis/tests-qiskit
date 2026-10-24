<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Math\Matrix\MatrixOperator;

class RYGate extends AbstractSingleQubitGate
{
    public const NAME = 'ry';

    private float $theta;

    public function __construct(float $theta, ?MatrixOperator $mo = null)
    {
        $this->theta = $theta;
        parent::__construct($mo);
    }

    public function name(): string
    {
        return self::NAME;
    }

    public function theta(): float
    {
        return $this->theta;
    }

    public function inverse(): self
    {
        return new self(-$this->theta, $this->mo);
    }

    protected function label(): string
    {
        return sprintf('%s(%.3g)', self::NAME, $this->theta);
    }

    protected function buildMatrix(): NDArray
    {
        $c = cos($this->theta / 2.0);
        $s = sin($this->theta / 2.0);
        return $this->la->array([[$c, -$s], [$s, $c]], dtype: NDArray::complex64);
    }
}
