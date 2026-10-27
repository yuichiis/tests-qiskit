<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;
use function Rindow\Math\Matrix\C;

class RZGate extends AbstractSingleQubitGate
{
    public const NAME = 'rz';

    private float $theta;

    public function __construct(float $theta, ?object $la=null)
    {
        $this->theta = $theta;
        parent::__construct(la:$la);
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
        return new self(-$this->theta, la:$this->la);
    }

    protected function label(): string
    {
        return sprintf('%s(%.3g)', self::NAME, $this->theta);
    }

    protected function buildMatrix(): NDArray
    {
        $c = cos($this->theta / 2.0);
        $s = sin($this->theta / 2.0);
        return $this->la->array(
            [[C($c, -$s), 0], [0, C($c, $s)]],
            dtype: NDArray::complex64
        );
    }
}
