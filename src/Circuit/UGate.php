<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Math\Matrix\MatrixOperator;
use function Rindow\Math\Matrix\C;

class UGate extends AbstractSingleQubitGate
{
    public const NAME = 'u';

    private float $theta;
    private float $phi;
    private float $lambda;

    public function __construct(float $theta, float $phi, float $lambda, ?MatrixOperator $mo = null)
    {
        $this->theta = $theta;
        $this->phi = $phi;
        $this->lambda = $lambda;
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

    public function phi(): float
    {
        return $this->phi;
    }

    public function lambda(): float
    {
        return $this->lambda;
    }

    public function inverse(): self
    {
        return new self(-$this->theta, -$this->lambda, -$this->phi, $this->mo);
    }

    protected function label(): string
    {
        return sprintf('%s(%.3g,%.3g,%.3g)', self::NAME, $this->theta, $this->phi, $this->lambda);
    }

    protected function buildMatrix(): NDArray
    {
        $c = cos($this->theta / 2.0);
        $s = sin($this->theta / 2.0);
        return $this->mo->array([
            [$c, C(-$s * cos($this->lambda), -$s * sin($this->lambda))],
            [
                C($s * cos($this->phi), $s * sin($this->phi)),
                C($c * cos($this->phi + $this->lambda), $c * sin($this->phi + $this->lambda)),
            ],
        ], dtype: NDArray::complex128);
    }
}
