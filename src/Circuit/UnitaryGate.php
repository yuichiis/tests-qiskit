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

    /** @var array<array{float,float}> Row-major complex [real, imaginary] pairs. */
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
        $out = array_fill(0, $size * $size, [0.0, 0.0]);
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
                    $out[$dst] = $this->matrix[$src];
                } elseif ($rowTarget === $colTarget) {
                    $out[$dst] = [1.0, 0.0];
                }
            }
        }
        return self::fromPairs($out);
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

    /** @return array<array{float,float}> Row-major complex [real, imaginary] pairs. */
    public function parameters(): array
    {
        return $this->matrix;
    }

    private static function fromPairs(array $pairs): self
    {
        $dim = (int)sqrt(count($pairs));
        $gate = new self(array_fill(0, $dim, array_fill(0, $dim, 0.0)));
        $gate->matrix = $pairs;
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
                if (is_array($value)) {
                    $out[] = [
                        (float)($value['re'] ?? $value[0] ?? 0.0),
                        (float)($value['im'] ?? $value[1] ?? 0.0),
                    ];
                } elseif (is_object($value) && method_exists($value, 'getReal')) {
                    $out[] = [(float)$value->getReal(), (float)$value->getImag()];
                } elseif (is_object($value) && (isset($value->real) || isset($value->re))) {
                    $out[] = [
                        (float)($value->real ?? $value->re ?? 0.0),
                        (float)($value->imag ?? $value->im ?? 0.0),
                    ];
                } else {
                    $out[] = [(float)$value, 0.0];
                }
            }
        }
        return $out;
    }

}
