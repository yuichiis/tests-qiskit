<?php
declare(strict_types=1);

namespace Rindow\Qubit\QuantumInfo;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Qubit\Math\LinearAlgebraProvider;
use Rindow\Qubit\Math\Precision;
use Rindow\Qubit\QuantumCircuit;
use function Rindow\Math\Matrix\C;

/**
 * qiskit.quantum_info.SparsePauliOp のPHP版.
 *
 * ```php
 * $h = SparsePauliOp::fromList([
 *     ['II', -1.052373245772859],
 *     ['ZZ', -0.01128010425623538],
 * ]);
 * $energy = $h->expectationValue($statevector);
 * ```
 */
class SparsePauliOp
{
    /** @var array<int,array{label:string,coeff:float}> */
    private array $terms = [];
    private int $numQubits = 0;
    /**
     * dtype ごとの行列キャッシュ. VQE のように期待値を繰り返し計算する
     * とき、項の展開と dtype への丸めを毎回やり直さない。
     *
     * @var array<int,NDArray>
     */
    private array $matrixCache = [];
    /**
     * dtype ごとの上三角パックド行列キャッシュ (BLAS hpmv 用).
     * n = 2^qubits として n(n+1)/2 要素だけ保持するため、密行列の
     * ほぼ半分のメモリで済む。
     *
     * @var array<int,NDArray>
     */
    private array $packedCache = [];

    /** @param array<int,array{0:string,1:float}> $terms */
    public function __construct(array $terms = [])
    {
        foreach ($terms as $term) {
            if (!is_array($term) || !array_key_exists(0, $term) || !array_key_exists(1, $term)) {
                throw new InvalidArgumentException('Each term must be [label, coefficient].');
            }
            $this->add((string)$term[0], (float)$term[1]);
        }
    }

    /** @param array<int,array{0:string,1:float}> $terms */
    public static function fromList(array $terms): self
    {
        return new self($terms);
    }

    public static function from_list(array $terms): self
    {
        return self::fromList($terms);
    }

    /** @param array<int,array{0:string,1:float}> $terms */
    public static function fromTerms(array $terms): self
    {
        return new self($terms);
    }

    private function add(string $label, float $coeff): void
    {
        $pauli = new Pauli($label);
        $numQubits = $pauli->numQubits();
        if ($this->numQubits === 0) {
            $this->numQubits = $numQubits;
        } elseif ($this->numQubits !== $numQubits) {
            throw new InvalidArgumentException(
                "Pauli label length mismatch: expected {$this->numQubits}, got {$numQubits}."
            );
        }
        if ($coeff == 0.0) {
            return;
        }
        $this->terms[] = ['label' => $label, 'coeff' => $coeff];
    }

    public function numQubits(): int
    {
        return $this->numQubits;
    }

    /** @deprecated use numQubits() */
    public function num_qubits(): int
    {
        return $this->numQubits;
    }

    public function numTerms(): int
    {
        return count($this->terms);
    }

    /**
     * @return array<int,array{label:string,coeff:float}>
     */
    public function terms(): array
    {
        return $this->terms;
    }

    /**
     * qiskit の to_list() 互換。
     *
     * @return array<int,array{0:string,1:float}>
     */
    public function toList(): array
    {
        $out = [];
        foreach ($this->terms as $term) {
            $out[] = [$term['label'], $term['coeff']];
        }
        return $out;
    }

    /** @return array<int,array{0:string,1:float}> */
    public function to_list(): array
    {
        return $this->toList();
    }

    /**
     * ハミルトニアン行列. 基底順は qubit 0 = LSB。
     *
     * dtype を省略すると Precision::complexDtype() を使う (既定 complex128)。
     *
     * @return NDArray shape [2^n, 2^n]
     */
    public function toMatrix(?int $dtype = null): NDArray
    {
        $dtype ??= Precision::complexDtype();
        // 呼び出し側が行列を書き換えてもキャッシュを壊さないようコピーを返す。
        return LinearAlgebraProvider::get()->copy($this->matrixForDtype($dtype));
    }

    /**
     * dtype ごとの行列を構築してキャッシュする. 内部利用専用.
     */
    private function matrixForDtype(int $dtype): NDArray
    {
        if ($dtype !== NDArray::complex64 && $dtype !== NDArray::complex128) {
            throw new InvalidArgumentException('dtype must be complex64 or complex128.');
        }
        if (isset($this->matrixCache[$dtype])) {
            return $this->matrixCache[$dtype];
        }
        if ($this->numQubits === 0) {
            throw new InvalidArgumentException('Empty SparsePauliOp has no matrix.');
        }
        $dim = 1 << $this->numQubits;
        $la = LinearAlgebraProvider::get();
        $matrix = $la->zeros([$dim, $dim], $dtype);
        foreach ($this->terms as $term) {
            $pauli = new Pauli($term['label']);
            // 複素配列では axpy の alpha も複素数でなければならない。
            // 負の係数はそのまま減算 (H = ... - c P) になる。
            $matrix = $la->axpy(
                $pauli->toMatrix($dtype),
                $matrix,
                C((float)$term['coeff'], 0.0),
            );
        }
        $this->matrixCache[$dtype] = $matrix;
        return $matrix;
    }

    /**
     * dtype ごとの上三角パックド行列を構築してキャッシュする.
     *
     * BLAS (hpmv) の RowMajor/Upper 形式に合わせ、行 i の j>=i を
     * 行順に並べる。H はエルミートなので下三角は保持しない。
     */
    private function packedMatrixForDtype(int $dtype): NDArray
    {
        if ($dtype !== NDArray::complex64 && $dtype !== NDArray::complex128) {
            throw new InvalidArgumentException('dtype must be complex64 or complex128.');
        }
        if (isset($this->packedCache[$dtype])) {
            return $this->packedCache[$dtype];
        }
        if ($this->numQubits === 0) {
            throw new InvalidArgumentException('Empty SparsePauliOp has no matrix.');
        }
        $dim = 1 << $this->numQubits;
        $la = LinearAlgebraProvider::get();
        $packed = $la->zeros([intdiv($dim * ($dim + 1), 2)], $dtype);
        foreach ($this->terms as $term) {
            $pauli = new Pauli($term['label']);
            $coeff = (float)$term['coeff'];
            for ($col = 0; $col < $dim; $col++) {
                [$row, $phase] = $pauli->applyToBasisState($col);
                // 上三角だけを格納する。パウリ行列はエルミートかつ
                // 対合 (P^2 = I) なので、下三角の成分は必ず別の列
                // (col = row) を処理するときに上三角側で現れる。
                if ($row > $col) {
                    continue;
                }
                $index = self::packedIndex($row, $col, $dim);
                $value = $packed[$index];
                $packed[$index] = C(
                    (float)($value->real ?? 0.0) + $coeff * $phase->real,
                    (float)($value->imag ?? 0.0) + $coeff * $phase->imag,
                );
            }
        }
        $this->packedCache[$dtype] = $packed;
        return $packed;
    }

    /** 上三角パックド行列 (RowMajor) の (row, col) の位置. row <= col. */
    private static function packedIndex(int $row, int $col, int $dim): int
    {
        return $row * $dim + $col - intdiv($row * ($row + 1), 2);
    }

    /** @return array<int,array<int,object{real:float,imag:float}>> */
    public function toArray(): array
    {
        return $this->toMatrix()->toArray();
    }

    /**
     * <ψ|H|ψ> を NDArray/BLAS で計算する.
     *
     * 状態ベクトルを PHP の float 配列に戻さず、エルミート・パックド
     * 行列ベクトル積 (hpmv) と複素内積 (dot, 共役あり) だけで期待値を
     * 求める。H は上三角だけ保持するので密行列よりメモリが少ない。
     *
     * @param Statevector|QuantumCircuit|NDArray|array<int,object{real:float,imag:float}> $state
     */
    public function expectationValue(Statevector|QuantumCircuit|NDArray|array $state): float
    {
        if ($state instanceof Statevector) {
            $ndarray = $state->ndarray();
        } elseif ($state instanceof QuantumCircuit) {
            $ndarray = (new Statevector($state))->ndarray();
        } elseif ($state instanceof NDArray) {
            $ndarray = $state;
        } else {
            $ndarray = (new Statevector($state))->ndarray();
        }
        return $this->expectationValueNDArray($ndarray);
    }

    /**
     * 状態ベクトル NDArray から <ψ|H|ψ> を計算する.
     */
    public function expectationValueNDArray(NDArray $state): float
    {
        $la = LinearAlgebraProvider::get();
        $dtype = $state->dtype();
        if ($dtype !== NDArray::complex64 && $dtype !== NDArray::complex128) {
            throw new InvalidArgumentException(
                'Statevector NDArray must have a complex dtype.'
            );
        }
        $shape = $state->shape();
        if (count($shape) !== 1 || $shape[0] !== (1 << $this->numQubits)) {
            throw new InvalidArgumentException(
                'Statevector qubit count does not match the operator.'
            );
        }

        // Hψ をエルミート・パックド行列積で作り、ψ†(Hψ) を BLAS の
        // 内積で求める。上三角だけを保持しているので gemv より省メモリ。
        $packed = $this->packedMatrixForDtype($dtype);
        $hpsi = $la->hpmv($packed, $state);
        $value = $la->scalar($la->dot($state, $hpsi));
        $real = (float)($value->real ?? 0.0);
        $imag = (float)($value->imag ?? 0.0);
        if (abs($imag) > 1.0e-5) {
            // エルミート演算子なら虚部はゼロになる。数値誤差だけを無視する。
            throw new InvalidArgumentException(
                "Expectation value has a non-negligible imaginary part: {$imag}"
            );
        }
        return $real;
    }

    /** @param Statevector|QuantumCircuit|NDArray|array<int,object{real:float,imag:float}> $state */
    public function expectation_value(Statevector|QuantumCircuit|NDArray|array $state): float
    {
        return $this->expectationValue($state);
    }
}
