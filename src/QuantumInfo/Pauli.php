<?php
declare(strict_types=1);

namespace Rindow\Qubit\QuantumInfo;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Qubit\Math\LinearAlgebraProvider;
use Rindow\Qubit\Math\Precision;
use function Rindow\Math\Matrix\C;

/**
 * qiskit.quantum_info.Pauli のPHP版.
 *
 * ラベルは qiskit と同じく左端が最上位qubit、右端が qubit 0。
 * 状態ベクトルの基底は qubit 0 = LSB。
 *
 * ```php
 * $p = new Pauli('IZ');          // qubit1 に I, qubit0 に Z
 * $m = $p->toMatrix();           // NDArray shape [4,4] (既定 complex128)
 * ```
 */
class Pauli
{
    private string $label;
    private int $numQubits;

    public function __construct(string $label)
    {
        $label = strtoupper($label);
        if ($label === '') {
            throw new InvalidArgumentException('Pauli label must not be empty.');
        }
        $valid = ['I' => true, 'X' => true, 'Y' => true, 'Z' => true];
        for ($i = 0, $len = strlen($label); $i < $len; $i++) {
            if (!isset($valid[$label[$i]])) {
                throw new InvalidArgumentException("Invalid Pauli character: {$label[$i]}");
            }
        }
        $this->label = $label;
        $this->numQubits = strlen($label);
    }

    public function label(): string
    {
        return $this->label;
    }

    public function numQubits(): int
    {
        return $this->numQubits;
    }

    /** qubit index $qubit (0 = LSB) のパウリ文字 */
    public function pauliAt(int $qubit): string
    {
        $this->checkQubit($qubit);
        return $this->label[$this->numQubits - 1 - $qubit];
    }

    /**
     * 計算基底 |index> に作用させたときの (遷移先index, 位相)。
     *
     * @return array{int,object{real:float,imag:float}}
     */
    public function applyToBasisState(int $index): array
    {
        $dim = 1 << $this->numQubits;
        if ($index < 0 || $index >= $dim) {
            throw new InvalidArgumentException("Basis index out of range: {$index}");
        }
        $target = $index;
        $phaseReal = 1.0;
        $phaseImag = 0.0;
        for ($k = 0; $k < $this->numQubits; $k++) {
            $q = $this->numQubits - 1 - $k;
            $char = $this->label[$k];
            $bit = ($index >> $q) & 1;
            switch ($char) {
                case 'I':
                    break;
                case 'X':
                    $target ^= 1 << $q;
                    break;
                case 'Y':
                    $target ^= 1 << $q;
                    // Y|0> = i|1>, Y|1> = -i|0>
                    [$phaseReal, $phaseImag] = self::multiply(
                        $phaseReal,
                        $phaseImag,
                        $bit === 0 ? 0.0 : -0.0,
                        $bit === 0 ? 1.0 : -1.0,
                    );
                    break;
                case 'Z':
                    if ($bit === 1) {
                        $phaseReal = -$phaseReal;
                        $phaseImag = -$phaseImag;
                    }
                    break;
            }
        }
        return [$target, (object)['real' => $phaseReal, 'imag' => $phaseImag]];
    }

    /**
     * パウリ行列を返す. 基底順は qubit 0 = LSB。
     *
     * dtype を省略すると Precision::complexDtype() を使う (既定 complex128)。
     *
     * @return NDArray shape [2^n, 2^n]
     */
    public function toMatrix(?int $dtype = null): NDArray
    {
        $dtype ??= Precision::complexDtype();
        if ($dtype !== NDArray::complex64 && $dtype !== NDArray::complex128) {
            throw new InvalidArgumentException('dtype must be complex64 or complex128.');
        }
        $la = LinearAlgebraProvider::get();
        $dim = 1 << $this->numQubits;
        $matrix = $la->zeros([$dim, $dim], $dtype);
        for ($col = 0; $col < $dim; $col++) {
            [$row, $phase] = $this->applyToBasisState($col);
            $matrix[$row][$col] = C($phase->real, $phase->imag);
        }
        return $matrix;
    }

    /** @return array<int,array<int,object{real:float,imag:float}>> */
    public function toArray(): array
    {
        return $this->toMatrix()->toArray();
    }

    public function __toString(): string
    {
        return $this->label;
    }

    private function checkQubit(int $qubit): void
    {
        if ($qubit < 0 || $qubit >= $this->numQubits) {
            throw new InvalidArgumentException("Qubit index out of range: {$qubit}");
        }
    }

    /** @return array{float,float} */
    private static function multiply(
        float $ar,
        float $ai,
        float $br,
        float $bi
    ): array {
        return [
            $ar * $br - $ai * $bi,
            $ar * $bi + $ai * $br,
        ];
    }
}
