<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Math\Matrix\MatrixOperator;
use function Rindow\Math\Matrix\C;

class PhaseGate extends AbstractSingleQubitGate
{
    public const NAME = 'p';

    private float $lambda;

    public function __construct(float $lambda, ?MatrixOperator $mo = null)
    {
        $this->lambda = $lambda;
        parent::__construct($mo);
    }

    public function name(): string
    {
        return self::NAME;
    }

    public function lambda(): float
    {
        return $this->lambda;
    }

    public function inverse(): self
    {
        return new self(-$this->lambda, $this->mo);
    }

    protected function label(): string
    {
        return sprintf('%s(%.3g)', self::NAME, $this->lambda);
    }

    protected function buildMatrix(): NDArray
    {
        return $this->la->array(
            [[1, 0], [0, C(cos($this->lambda), sin($this->lambda))]],
            dtype: NDArray::complex64
        );
    }
}
