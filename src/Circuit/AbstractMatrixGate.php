<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Qubit\CircuitGate;
use Rindow\Qubit\Math\LinearAlgebraProvider;
use Rindow\Qubit\Math\Precision;

abstract class AbstractMatrixGate implements CircuitGate
{
    /** @var NDArray shape [dim, dim], row-major. dtype は Precision::complexDtype()。 */
    protected NDArray $matrix;
    protected object $la;

    public function __construct(?object $la=null)
    {
        $this->la = $la ?? LinearAlgebraProvider::get();
        $this->matrix = $this->buildMatrix();
        $dtype = $this->complexDtype();
        if ($this->matrix->dtype() !== $dtype) {
            // matlib の astype は complex64<->complex128 に未対応のため
            // PHP 配列経由で変換する。
            $this->matrix = $this->la->array(
                $this->matrix->toArray(),
                dtype: $dtype,
            );
        }
    }

    /** ゲート行列の複素 dtype (Precision の設定に従う). */
    protected function complexDtype(): int
    {
        return Precision::complexDtype();
    }

    abstract protected function buildMatrix(): NDArray;

    public function dimension(): int
    {
        return $this->matrix->shape()[0];
    }

    /** @return NDArray shape [dim, dim], dtype は Precision::complexDtype() */
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
