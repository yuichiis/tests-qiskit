<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Math\Matrix\Complex;
use Rindow\Qubit\CircuitGate;
use Rindow\Qubit\Math\LinearAlgebraProvider;
use function Rindow\Math\Matrix\C;

class UnitaryGate implements CircuitGate
{
    public const NAME = 'unitary';

    /** @var NDArray complex64, shape [dim, dim], row-major. */
    private NDArray $matrix;
    private object $la;

    public function __construct(
        array|NDArray $matrix,
        ?object $la=null
    ) {
        $this->la = $la ?? LinearAlgebraProvider::get();
        $this->matrix = $this->normalizeMatrix($matrix);
    }

    public function name(): string
    {
        return self::NAME;
    }

    public function draw(array $qubits): array
    {
        $cells = [];
        $symbol = count($qubits) === 1 ? self::NAME : 'U';
        foreach ($qubits as $q) {
            $cells[$q] = $symbol;
        }
        $span = range(min($qubits), max($qubits));
        for ($q = min($qubits) + 1; $q < max($qubits); $q++) {
            $cells[$q] = '│';
        }
        return ['cells' => $cells, 'span' => $span];
    }

    public function control(int $numControls = 1): self
    {
        if ($numControls < 1) {
            throw new InvalidArgumentException('numControls must be at least 1.');
        }
        $dim = $this->dimension();
        $controlDim = 1 << $numControls;
        $size = $dim * $controlDim;
        $out = $this->la->zeros([$size * $size], NDArray::complex64);
        $src = $this->matrix->reshape([$dim * $dim]);
        for ($row = 0; $row < $size; $row++) {
            for ($col = 0; $col < $size; $col++) {
                $rowControl = $row % $controlDim;
                $colControl = $col % $controlDim;
                if ($rowControl !== $colControl) continue;
                $rowTarget = intdiv($row, $controlDim);
                $colTarget = intdiv($col, $controlDim);
                $dst = $row * $size + $col;
                if ($rowControl === $controlDim - 1) {
                    $srcIndex = $rowTarget * $dim + $colTarget;
                    $out[$dst] = $src[$srcIndex];
                } elseif ($rowTarget === $colTarget) {
                    $out[$dst] = C(1.0, 0.0);
                }
            }
        }
        return new self($out->reshape([$size, $size]), $this->la);
    }

    public static function fredkin(): self
    {
        $la = LinearAlgebraProvider::get();
        $dim = 8;
        $matrix = $la->zeros([$dim * $dim], NDArray::complex64);
        for ($input = 0; $input < $dim; $input++) {
            $control = $input & 1;
            $a = ($input >> 1) & 1;
            $b = ($input >> 2) & 1;
            $output = $control ? ($control | ($b << 1) | ($a << 2)) : $input;
            $matrix[$output * $dim + $input] = C(1.0, 0.0);
        }
        return new self($matrix->reshape([$dim, $dim]), la:$la);
    }

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
        $flat = $this->matrix->reshape([$dim * $dim]);
        $out = [];
        for ($i = 0; $i < $dim * $dim; $i++) {
            $c = $flat[$i];
            $out[] = (object)['real' => (float)$c->real, 'imag' => (float)$c->imag];
        }
        return $out;
    }

    private function normalizeMatrix(array|NDArray $matrix): NDArray
    {
        if ($matrix instanceof NDArray) {
            $shape = $matrix->shape();
            if (count($shape) !== 2 || $shape[0] !== $shape[1]) {
                throw new InvalidArgumentException('unitary() requires an NxN matrix (N=2^k).');
            }
            $this->checkDimension($shape[0]);
            if ($matrix->dtype() !== NDArray::complex64) {
                return $this->la->array($matrix->toArray(), dtype: NDArray::complex64);
            }
            return $matrix;
        }

        $rows = array_values($matrix);
        $dim = count($rows);
        $this->checkDimension($dim);
        $flat = $this->la->zeros([$dim * $dim], NDArray::complex64);
        $idx = 0;
        foreach ($rows as $row) {
            if (!is_array($row) || count($row) !== $dim) {
                throw new InvalidArgumentException('unitary() requires an NxN matrix (N=2^k).');
            }
            foreach ($row as $value) {
                $flat[$idx++] = self::toComplex($value);
            }
        }
        return $flat->reshape([$dim, $dim]);
    }

    private function checkDimension(int $dim): void
    {
        if ($dim < 2 || ($dim & ($dim - 1)) !== 0) {
            throw new InvalidArgumentException('unitary() requires an NxN matrix (N=2^k).');
        }
    }

    private static function toComplex(mixed $value): Complex
    {
        if (is_object($value)) {
            if (method_exists($value, 'getReal')) {
                return C((float)$value->getReal(), (float)$value->getImag());
            }
            if (property_exists($value, 'real') || isset($value->real)) {
                return C((float)$value->real, (float)($value->imag ?? 0.0));
            }
            throw new InvalidArgumentException(
                'unitary() elements must be real numbers or complex objects with real/imag.'
            );
        }
        if (is_array($value)) {
            if (array_key_exists('real', $value) || array_key_exists('imag', $value)) {
                return C(
                    (float)($value['real'] ?? 0.0),
                    (float)($value['imag'] ?? 0.0),
                );
            }
            throw new InvalidArgumentException(
                'unitary() does not accept plain [re, im] pairs; use real/imag.'
            );
        }
        return C((float)$value, 0.0);
    }
}
