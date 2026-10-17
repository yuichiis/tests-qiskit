<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use Rindow\Qubit\CircuitUtils;

class UnitaryGate extends AbstractGate
{
    use CircuitUtils;

    public const NAME = 'unitary';

    /**
     * @param array<int,int> $qubits
     */
    public function __construct(
        private array $qubits,
    ) {
        parent::__construct();
    }

}