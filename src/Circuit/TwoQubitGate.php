<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

abstract class TwoQubitGate extends MatrixGate
{
    public function dimension(): int
    {
        return 4;
    }

    protected function drawTwoQubit(array $qubits, string $first, string $second): array
    {
        [$a, $b] = $qubits;
        $cells = [$a => $first, $b => $second];
        $span = range(min($a, $b), max($a, $b));
        for ($q = min($a, $b) + 1; $q < max($a, $b); $q++) {
            $cells[$q] = '│';
        }
        return ['cells' => $cells, 'span' => $span];
    }
}
