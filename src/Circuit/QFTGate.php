<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use function Rindow\Math\Matrix\C;

class QFTGate extends AbstractMatrixGate
{
    public const NAME = 'qft';
    public const NAME_INVERSE = 'qft_dg';

    private int $numQubits;
    private bool $inverse = false;

    public function __construct(
        int $numQubits,
        ?object $la=null,
    ) {
        if ($numQubits < 1 || $numQubits > 30) {
            throw new InvalidArgumentException('numQubits must be in [1, 30].');
        }
        $this->numQubits = $numQubits;
        parent::__construct($la);
    }

    public function inverse(): self
    {
        $inverse = clone $this;
        $inverse->inverse = !$this->inverse;
        $inverse->matrix = $inverse->buildMatrix();
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

    protected function buildMatrix(): NDArray
    {
        $dim = 1 << $this->numQubits;
        $scale = 1.0 / sqrt((float)$dim);
        $sign = $this->inverse ? -1.0 : 1.0;
        $flat = $this->la->zeros([$dim * $dim], NDArray::complex64);
        for ($j = 0; $j < $dim; $j++) {
            for ($k = 0; $k < $dim; $k++) {
                $angle = $sign * 2.0 * M_PI * $j * $k / $dim;
                $flat[$j * $dim + $k] = C($scale * cos($angle), $scale * sin($angle));
            }
        }
        return $flat->reshape([$dim, $dim]);
    }
}
