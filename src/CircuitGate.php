<?php
declare(strict_types=1);

namespace Rindow\Qubit;

interface CircuitGate
{
    /**
     * @param array<int> $qubits
     * @return array{cells:array<int,string>,span:array<int>}
     */
    public function draw(array $qubits): array;
}
