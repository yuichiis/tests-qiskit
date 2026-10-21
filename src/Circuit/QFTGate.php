<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Qubit\CircuitGate;
use Rindow\Math\Matrix\MatrixOperator;
use Rindow\Qubit\Math\MatrixOperatorProvider;
use function Rindow\Math\Matrix\C;

class QFTGate implements CircuitGate
{
    public const NAME = 'qft';
    public const NAME_INVERSE = 'qft_dg';

    private int $numQubits;
    private bool $inverse = false;
    /** @var NDArray complex128, shape [dim, dim], row-major. */
    private NDArray $matrix;
    private MatrixOperator $mo;

    public function __construct(int $numQubits, ?MatrixOperator $mo = null)
    {
        if ($numQubits < 1 || $numQubits > 30) {
            throw new InvalidArgumentException('numQubits must be in [1, 30].');
        }
        $this->mo = $mo ?? MatrixOperatorProvider::get();
        $this->numQubits = $numQubits;
        $this->matrix = $this->buildMatrix(false);
    }

    public function inverse(): self
    {
        $inverse = clone $this;
        $inverse->inverse = !$this->inverse;
        $inverse->matrix = $inverse->buildMatrix($inverse->inverse);
        return $inverse;
    }

    public function name(): string
    {
        return $this->inverse ? self::NAME_INVERSE : self::NAME;
    }

    public function draw(array $qubits): array
    {
        $symbol = $this->inverse ? 'QFT†' : 'QFT';
        $cells = [];
        foreach ($qubits as $q) {
            $cells[$q] = $symbol;
        }
        $span = range(min($qubits), max($qubits));
        for ($q = min($qubits) + 1; $q < max($qubits); $q++) {
            $cells[$q] = '│';
        }
        return ['cells' => $cells, 'span' => $span];
    }

    public function isInverse(): bool
    {
        return $this->inverse;
    }

    public function numQubits(): int
    {
        return $this->numQubits;
    }

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

    private function buildMatrix(bool $inverse): NDArray
    {
        $dim = 1 << $this->numQubits;
        $scale = 1.0 / sqrt((float)$dim);
        $sign = $inverse ? -1.0 : 1.0;
        $flat = $this->mo->zeros([$dim * $dim], NDArray::complex128);
        for ($j = 0; $j < $dim; $j++) {
            for ($k = 0; $k < $dim; $k++) {
                $angle = $sign * 2.0 * M_PI * $j * $k / $dim;
                $flat[$j * $dim + $k] = C($scale * cos($angle), $scale * sin($angle));
            }
        }
        return $flat->reshape([$dim, $dim]);
    }
}
