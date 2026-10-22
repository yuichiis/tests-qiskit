<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

abstract class AbstractSingleQubitGate extends AbstractMatrixGate
{
    public function dimension(): int
    {
        return 2;
    }

    protected function label(): string
    {
        return $this->name();
    }

    public function draw(array $qubits): array
    {
        $q = $qubits[0];
        return ['cells' => [$q => $this->label()], 'span' => [$q]];
    }
}
