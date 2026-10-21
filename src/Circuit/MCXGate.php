<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Qubit\CircuitGate;
use Rindow\Math\Matrix\MatrixOperator;
use Rindow\Qubit\Math\MatrixOperatorProvider;
use function Rindow\Math\Matrix\C;

class MCXGate implements CircuitGate
{
    public const NAME = 'mcx';

    private int $numCtrlQubits;
    private int $ctrlState;
    /** @var NDArray complex128, shape [dim, dim], row-major. */
    private NDArray $matrix;
    private MatrixOperator $mo;

    public function __construct(
        int $numCtrlQubits,
        int|string|null $ctrlState = null,
        ?MatrixOperator $mo = null
    ) {
        if ($numCtrlQubits < 1 || $numCtrlQubits > 29) {
            throw new InvalidArgumentException('numCtrlQubits must be in [1, 29].');
        }
        $this->mo = $mo ?? MatrixOperatorProvider::get();
        $this->numCtrlQubits = $numCtrlQubits;
        $this->ctrlState = self::normalizeCtrlState($ctrlState, $numCtrlQubits);
        $this->matrix = $this->buildMatrix();
    }

    public function name(): string
    {
        return self::NAME;
    }

    public function draw(array $qubits): array
    {
        $last = count($qubits) - 1;
        $cells = [];
        foreach ($qubits as $i => $q) {
            $cells[$q] = $i === $last ? '⊕' : '●';
        }
        $span = range(min($qubits), max($qubits));
        foreach ($span as $q) {
            if (!array_key_exists($q, $cells)) {
                $cells[$q] = '│';
            }
        }
        return ['cells' => $cells, 'span' => $span];
    }

    public function numCtrlQubits(): int
    {
        return $this->numCtrlQubits;
    }

    public function numQubits(): int
    {
        return $this->numCtrlQubits + 1;
    }

    public function ctrlState(): int
    {
        return $this->ctrlState;
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

    private function buildMatrix(): NDArray
    {
        $dim = 1 << $this->numQubits();
        $targetBit = 1 << $this->numCtrlQubits;
        $controlMask = $targetBit - 1;
        $flat = $this->mo->zeros([$dim * $dim], NDArray::complex128);
        for ($input = 0; $input < $dim; $input++) {
            $output = $input;
            if (($input & $controlMask) === $this->ctrlState) {
                $output = $input ^ $targetBit;
            }
            $flat[$output * $dim + $input] = C(1.0, 0.0);
        }
        return $flat->reshape([$dim, $dim]);
    }

    private static function normalizeCtrlState(int|string|null $ctrlState, int $numCtrlQubits): int
    {
        $limit = 1 << $numCtrlQubits;
        if ($ctrlState === null) {
            return $limit - 1;
        }
        if (is_string($ctrlState)) {
            if (strlen($ctrlState) !== $numCtrlQubits || preg_match('/^[01]+$/', $ctrlState) !== 1) {
                throw new InvalidArgumentException(
                    "ctrlState must be a {$numCtrlQubits}-bit binary string."
                );
            }
            return (int)bindec($ctrlState);
        }
        if ($ctrlState < 0 || $ctrlState >= $limit) {
            throw new InvalidArgumentException("ctrlState must be in [0, {$limit}).");
        }
        return $ctrlState;
    }
}
