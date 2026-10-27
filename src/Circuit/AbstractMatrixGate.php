<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Qubit\CircuitGate;
use Rindow\Qubit\Math\LinearAlgebraProvider;

abstract class AbstractMatrixGate implements CircuitGate
{
    /** @var NDArray complex64, shape [dim, dim], row-major. */
    protected NDArray $matrix;
    protected object $la;

    public function __construct(?object $la=null)
    {
        $this->la = $la ?? LinearAlgebraProvider::get();
        $this->matrix = $this->buildMatrix();
    }

    abstract protected function buildMatrix(): NDArray;

    public function dimension(): int
    {
        return $this->matrix->shape()[0];
    }

    /** @return NDArray complex64, shape [dim, dim] */
    public function matrix(): NDArray
    {
        return $this->matrix;
    }

    /** @return array<object{real:float,imag:float}> Row-major complex numbers. */
    public function parameters(): array
    {
        $dim = $this->dimension();
        $values = $this->matrix->reshape([$dim * $dim])->toArray();
        $out = [];
        foreach ($values as $c) {
            $out[] = (object)['real' => (float)$c->real, 'imag' => (float)$c->imag];
        }
        return $out;
    }
}
