<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Qubit\CircuitGate;
use Rindow\Math\Matrix\MatrixOperator;
use Rindow\Qubit\Math\MatrixOperatorProvider;

abstract class MatrixGate implements CircuitGate
{
    /** @var NDArray complex128, shape [dim, dim], row-major. */
    protected NDArray $matrix;
    protected MatrixOperator $mo;

    public function __construct(?MatrixOperator $mo = null)
    {
        $this->mo = $mo ?? MatrixOperatorProvider::get();
        $this->matrix = $this->buildMatrix();
    }

    abstract protected function buildMatrix(): NDArray;

    public function dimension(): int
    {
        return $this->matrix->shape()[0];
    }

    /** @return NDArray complex128, shape [dim, dim] */
    public function matrix(): NDArray
    {
        return $this->matrix;
    }

    /** @return array<object{real:float,imag:float}> Row-major complex numbers. */
    public function parameters(): array
    {
        $dim = $this->dimension();
        $flat = $this->matrix->reshape([$dim * $dim]);
        $out = [];
        for ($i = 0; $i < $dim * $dim; $i++) {
            $c = $flat[$i];
            $out[] = (object)['real' => (float)$c->real, 'imag' => (float)$c->imag];
        }
        return $out;
    }
}
