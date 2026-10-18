<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Rindow\Qubit\CircuitUtils;
use InvalidArgumentException;

use Rindow\Qubit\CircuitGate;

class UnitaryGate implements CircuitGate
{
    use CircuitUtils;

    public const NAME = 'unitary';

    /** @var array<object{real:float,imag:float}> Row-major complex numbers. */
    private array $matrix;

    public function __construct(array $matrix)
    {
        $this->matrix = self::matrixToComplexPairs($matrix);
    }

    public function control(int $numControls = 1): self
    {
        if ($numControls < 1) {
            throw new InvalidArgumentException('numControls must be at least 1.');
        }
        $dim = (int)sqrt(count($this->matrix));
        $controlDim = 1 << $numControls;
        $size = $dim * $controlDim;
        $out = [];
        for ($i = 0, $total = $size * $size; $i < $total; $i++) {
            $out[] = (object)['real' => 0.0, 'imag' => 0.0];
        }
        for ($row = 0; $row < $size; $row++) {
            for ($col = 0; $col < $size; $col++) {
                $rowControl = $row % $controlDim;
                $colControl = $col % $controlDim;
                if ($rowControl !== $colControl) continue;
                $rowTarget = intdiv($row, $controlDim);
                $colTarget = intdiv($col, $controlDim);
                $dst = $row * $size + $col;
                if ($rowControl === $controlDim - 1) {
                    $src = $rowTarget * $dim + $colTarget;
                    $out[$dst] = clone $this->matrix[$src];
                } elseif ($rowTarget === $colTarget) {
                    $out[$dst] = (object)['real' => 1.0, 'imag' => 0.0];
                }
            }
        }
        return self::fromComplex($out);
    }

    public static function fredkin(): self
    {
        $matrix = array_fill(0, 8, array_fill(0, 8, 0.0));
        for ($input = 0; $input < 8; $input++) {
            $control = $input & 1;
            $a = ($input >> 1) & 1;
            $b = ($input >> 2) & 1;
            $output = $control ? ($control | ($b << 1) | ($a << 2)) : $input;
            $matrix[$output][$input] = 1.0;
        }
        return new self($matrix);
    }

    /** @return array<object{real:float,imag:float}> Row-major complex numbers. */
    public function parameters(): array
    {
        return $this->matrix;
    }

    private static function fromComplex(array $values): self
    {
        $dim = (int)sqrt(count($values));
        $gate = new self(array_fill(0, $dim, array_fill(0, $dim, 0.0)));
        $gate->matrix = $values;
        return $gate;
    }

    private static function matrixToComplexPairs(array $matrix): array
    {
        $rows = array_values($matrix);
        $dim = count($rows);
        if ($dim < 2 || ($dim & ($dim - 1)) !== 0) {
            throw new InvalidArgumentException('unitary() requires an NxN matrix (N=2^k).');
        }
        $out = [];
        foreach ($rows as $row) {
            if (!is_array($row) || count($row) !== $dim) {
                throw new InvalidArgumentException('unitary() requires an NxN matrix (N=2^k).');
            }
            foreach ($row as $value) {
                $out[] = self::toComplex($value);
            }
        }
        return $out;
    }

    private static function toComplex(mixed $value): object
    {
        if (is_object($value)) {
            if (method_exists($value, 'getReal')) {
                return (object)[
                    'real' => (float)$value->getReal(),
                    'imag' => (float)$value->getImag(),
                ];
            }
            if (property_exists($value, 'real') || isset($value->real)) {
                return (object)[
                    'real' => (float)$value->real,
                    'imag' => (float)($value->imag ?? 0.0),
                ];
            }
            throw new InvalidArgumentException(
                'unitary() elements must be real numbers or complex objects with real/imag.'
            );
        }
        if (is_array($value)) {
            if (array_key_exists('real', $value) || array_key_exists('imag', $value)) {
                return (object)[
                    'real' => (float)($value['real'] ?? 0.0),
                    'imag' => (float)($value['imag'] ?? 0.0),
                ];
            }
            throw new InvalidArgumentException(
                'unitary() does not accept plain [re, im] pairs; use real/imag.'
            );
        }
        return (object)['real' => (float)$value, 'imag' => 0.0];
    }

}
