<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Qubit\CircuitGate;
use Rindow\Qubit\Math\LinearAlgebraProvider;
use function Rindow\Math\Matrix\C;

class MCPhaseGate implements CircuitGate
{
    public const NAME = 'mcphase';

    private float $lambda;
    private int $numCtrlQubits;
    private int $ctrlState;
    /** @var NDArray complex64, shape [dim, dim], row-major. */
    private NDArray $matrix;
    private object $la;

    public function __construct(
        float $lambda,
        int $numCtrlQubits,
        int|string|null $ctrlState = null,
        ?object $la=null
    ) {
        if ($numCtrlQubits < 1 || $numCtrlQubits > 29) {
            throw new InvalidArgumentException('numCtrlQubits must be in [1, 29].');
        }
        $this->la = $la ?? LinearAlgebraProvider::get();
        $this->lambda = $lambda;
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
            $cells[$q] = $i === $last ? 'P' : '●';
        }
        $span = range(min($qubits), max($qubits));
        foreach ($span as $q) {
            if (!array_key_exists($q, $cells)) {
                $cells[$q] = '│';
            }
        }
        return ['cells' => $cells, 'span' => $span];
    }

    public function inverse(): self
    {
        return new self(-$this->lambda, $this->numCtrlQubits, $this->ctrlState, $this->la);
    }

    public function lambda(): float
    {
        return $this->lambda;
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

    private function buildMatrix(): NDArray
    {
        $dim = 1 << $this->numQubits();
        $targetBit = 1 << $this->numCtrlQubits;
        $controlMask = $targetBit - 1;
        $phase = C(cos($this->lambda), sin($this->lambda));
        $flat = $this->la->zeros([$dim * $dim], NDArray::complex64);
        for ($input = 0; $input < $dim; $input++) {
            $value = C(1.0, 0.0);
            if (($input & $controlMask) === $this->ctrlState && ($input & $targetBit) !== 0) {
                $value = $phase;
            }
            $flat[$input * $dim + $input] = $value;
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
