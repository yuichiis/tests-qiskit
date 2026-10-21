<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Qubit\CircuitGate;
use Rindow\Math\Matrix\MatrixOperator;
use Rindow\Qubit\Math\MatrixOperatorProvider;

class HGate implements CircuitGate
{
    public const NAME = 'h';

    /** @var NDArray complex128, shape [2, 2], row-major. */
    private NDArray $matrix;
    private MatrixOperator $mo;

    public function __construct(?MatrixOperator $mo = null)
    {
        $this->mo = $mo ?? MatrixOperatorProvider::get();
        $s = 1.0 / sqrt(2.0);
        $this->matrix = $this->mo->array(
            [[$s, $s], [$s, -$s]],
            dtype: NDArray::complex128
        );
    }

    public function name(): string
    {
        return self::NAME;
    }

    public function draw(array $qubits): array
    {
        $q = $qubits[0];
        return ['cells' => [$q => self::NAME], 'span' => [$q]];
    }

    public function inverse(): self
    {
        return new self($this->mo);
    }

    public function dimension(): int
    {
        return 2;
    }

    /** @return NDArray complex128, shape [2, 2] */
    public function matrix(): NDArray
    {
        return $this->matrix;
    }

    /** @return array<object{real:float,imag:float}> Row-major complex numbers. */
    public function parameters(): array
    {
        $flat = $this->matrix->reshape([4]);
        $out = [];
        for ($i = 0; $i < 4; $i++) {
            $c = $flat[$i];
            $out[] = (object)['real' => (float)$c->real, 'imag' => (float)$c->imag];
        }
        return $out;
    }
}
