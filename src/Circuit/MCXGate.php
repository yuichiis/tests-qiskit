<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use function Rindow\Math\Matrix\C;

class MCXGate extends AbstractMatrixGate
{
    public const NAME = 'mcx';

    private int $numCtrlQubits;
    private int $ctrlState;

    public function __construct(
        int $numCtrlQubits,
        int|string|null $ctrlState = null,
        ?object $la=null
    ) {
        if ($numCtrlQubits < 1 || $numCtrlQubits > 29) {
            throw new InvalidArgumentException('numCtrlQubits must be in [1, 29].');
        }
        $this->numCtrlQubits = $numCtrlQubits;
        $this->ctrlState = self::normalizeCtrlState($ctrlState, $numCtrlQubits);
        parent::__construct($la);
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

    protected function buildMatrix(): NDArray
    {
        $dim = 1 << $this->numQubits();
        $targetBit = 1 << $this->numCtrlQubits;
        $controlMask = $targetBit - 1;
        $flat = $this->la->zeros([$dim * $dim], NDArray::complex64);
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
