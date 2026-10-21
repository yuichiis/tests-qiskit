<?php
declare(strict_types=1);

namespace Rindow\Qubit;

use Interop\Polite\Math\Matrix\NDArray;

interface CircuitGate
{
    public function name(): string;

    public function dimension(): int;

    /** @return NDArray complex128, shape [dim, dim] */
    public function matrix(): NDArray;

    /** @return array<object{real:float,imag:float}> Row-major complex numbers. */
    public function parameters(): array;

    /**
     * @param array<int> $qubits
     * @return array{cells:array<int,string>,span:array<int>}
     */
    public function draw(array $qubits): array;
}
