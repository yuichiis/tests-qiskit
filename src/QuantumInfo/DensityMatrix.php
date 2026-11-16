<?php
declare(strict_types=1);

namespace Rindow\Qubit\QuantumInfo;

use Countable;
use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Math\Matrix\Complex;
use Rindow\Qubit\Math\IntMath;
use Rindow\Qubit\Math\LinearAlgebraProvider;
use Rindow\Qubit\Math\Precision;
use Rindow\Qubit\Quantum\StatevectorSimulator;
use Rindow\Qubit\QuantumCircuit;
use function Rindow\Math\Matrix\C;

/**
 * qiskit.quantum_info.DensityMatrix のPHP版.
 *
 * 内部では complex の NDArray (shape [2^n, 2^n]) を保持し、
 * 状態ベクトルからの生成には BLAS の gemm (A A^H) を使う。
 * 混合状態は fromEnsemble() で gemm の beta 累積を使って組み立てる。
 *
 * ```php
 * use Rindow\Qubit\QuantumInfo\DensityMatrix;
 *
 * $rho = DensityMatrix::fromStatevector($state);
 * $reduced = $rho->partialTrace([1]);   // qubit 1 をトレースアウト
 * $purity = $rho->purity();
 * ```
 */
class DensityMatrix implements Countable
{
    private NDArray $ndarray;
    private int $numQubits;
    private int $dim;

    /**
     * @param NDArray|array|Statevector|StatevectorSimulator|QuantumCircuit|DensityMatrix $data
     */
    public function __construct(
        NDArray|array|Statevector|StatevectorSimulator|QuantumCircuit|self $data
    ) {
        if ($data instanceof self) {
            $this->assignMatrix($data->ndarray);
            return;
        }
        if ($data instanceof Statevector
            || $data instanceof StatevectorSimulator
            || $data instanceof QuantumCircuit
        ) {
            $this->assignMatrix(self::densityFromVector(self::vectorNDArray($data)));
            return;
        }
        if ($data instanceof NDArray) {
            $shape = $data->shape();
            if (count($shape) === 1) {
                $this->assignMatrix(self::densityFromVector($data));
                return;
            }
            if (count($shape) === 2 && $shape[1] === 1) {
                $vector = $data->reshape([$shape[0]]);
                $this->assignMatrix(self::densityFromVector($vector));
                return;
            }
            $this->assignMatrix($data);
            return;
        }

        if ($data === []) {
            throw new InvalidArgumentException('The density matrix must not be empty.');
        }
        $first = $data[array_key_first($data)];
        if (is_array($first) || $first instanceof \ArrayObject) {
            $la = LinearAlgebraProvider::get();
            $this->assignMatrix($la->array($data, dtype: Precision::complexDtype()));
            return;
        }
        // 1次元の配列は状態ベクトルとして扱い、|ψ><ψ| を作る。
        $this->assignMatrix(self::densityFromVector(self::vectorNDArray($data)));
    }

    /**
     * 状態ベクトルから密度行列を作る.
     *
     * $qargs を省略すると |ψ><ψ| 全体を返す。$qargs を指定した場合は
     * その qubit だけを残した縮約密度行列を、ψ.reshape + gemm で直接作る。
     * $qargs の先頭が結果の最下位ビットになる (probabilities() と同じ順序)。
     *
     * @param Statevector|StatevectorSimulator|QuantumCircuit|NDArray|array $state
     * @param array<int>|null $qargs 残す qubit のリスト
     */
    public static function fromStatevector(
        Statevector|StatevectorSimulator|QuantumCircuit|NDArray|array $state,
        ?array $qargs = null
    ): self {
        $vector = self::vectorNDArray($state);
        if ($qargs === null) {
            return new self(self::densityFromVector($vector));
        }
        return new self(self::reducedDensityFromVector($vector, $qargs));
    }

    /**
     * 状態ベクトルのアンサンブルから混合状態を作る.
     *
     * ρ = (1/N) Σ_k |ψ_k><ψ_k| を gemm の beta 累積で計算する。
     * $normalize が false のときは和のまま返す。
     *
     * @param array<int,Statevector|StatevectorSimulator|QuantumCircuit|NDArray|array> $states
     */
    public static function fromEnsemble(array $states, bool $normalize = true): self
    {
        if ($states === []) {
            throw new InvalidArgumentException('The ensemble must not be empty.');
        }
        $la = LinearAlgebraProvider::get();
        $rho = null;
        $expectedDim = null;
        $dtype = null;
        foreach ($states as $state) {
            $vector = self::vectorNDArray($state);
            $dim = $vector->shape()[0];
            if ($expectedDim === null) {
                $expectedDim = $dim;
            } elseif ($expectedDim !== $dim) {
                throw new InvalidArgumentException(
                    'All statevectors in the ensemble must have the same dimension.'
                );
            }
            if ($dtype === null) {
                $dtype = $vector->dtype();
            } elseif ($vector->dtype() !== $dtype) {
                // 精度の異なる状態ベクトルを混ぜると計算結果が壊れるため、
                // 昇格や丸めはせずエラーにする。
                throw new InvalidArgumentException(
                    'All statevectors in the ensemble must have the same precision.'
                );
            }
            $column = $vector->reshape([$dim, 1]);
            $rho = self::projector($column, $rho);
        }
        /** @var NDArray $rho */
        if ($normalize) {
            $la->scal(C(1.0 / count($states), 0.0), $rho);
        }
        return new self($rho);
    }

    /**
     * |i><i| を作る.
     */
    public static function fromInt(int $index, int $numQubits): self
    {
        if ($numQubits < 0) {
            throw new InvalidArgumentException('The number of qubits must not be negative.');
        }
        $dim = 1 << $numQubits;
        if ($index < 0 || $index >= $dim) {
            throw new InvalidArgumentException("Basis index out of range: {$index}");
        }
        $la = LinearAlgebraProvider::get();
        $dtype = Precision::complexDtype();
        $matrix = $la->zeros([$dim, $dim], $dtype);
        $matrix[$index][$index] = C(1.0, 0.0);
        return new self($matrix);
    }

    private function assignMatrix(NDArray $matrix): void
    {
        $la = LinearAlgebraProvider::get();
        $shape = $matrix->shape();
        if (count($shape) !== 2 || $shape[0] !== $shape[1]) {
            throw new InvalidArgumentException(
                'The density matrix must be a square matrix.'
            );
        }
        $dim = $shape[0];
        self::assertDimension($dim);
        $dtype = $matrix->dtype();
        if ($dtype !== NDArray::complex64 && $dtype !== NDArray::complex128) {
            $matrix = $la->array($matrix->toArray(), dtype: Precision::complexDtype());
        }
        $this->ndarray = $la->copy($matrix);
        $this->dim = $dim;
        $this->numQubits = $dim === 1 ? 0 : (int)round(log($dim, 2));
    }

    private static function assertDimension(int $dim): void
    {
        if ($dim < 1 || ($dim & ($dim - 1)) !== 0) {
            throw new InvalidArgumentException(
                'The dimension of the density matrix must be a power of two.'
            );
        }
    }

    /**
     * @param Statevector|StatevectorSimulator|QuantumCircuit|NDArray|array $state
     */
    private static function vectorNDArray(
        Statevector|StatevectorSimulator|QuantumCircuit|NDArray|array $state
    ): NDArray {
        $la = LinearAlgebraProvider::get();
        if ($state instanceof Statevector) {
            return $state->ndarray();
        }
        if ($state instanceof StatevectorSimulator) {
            return (new Statevector($state))->ndarray();
        }
        if ($state instanceof QuantumCircuit) {
            return (new Statevector($state))->ndarray();
        }
        if ($state instanceof NDArray) {
            $shape = $state->shape();
            if (count($shape) !== 1) {
                throw new InvalidArgumentException(
                    'The statevector NDArray must be one-dimensional.'
                );
            }
            self::assertDimension($shape[0]);
            $dtype = $state->dtype();
            if ($dtype !== NDArray::complex64 && $dtype !== NDArray::complex128) {
                return $la->array($state->toArray(), dtype: Precision::complexDtype());
            }
            return $la->copy($state);
        }
        if ($state === []) {
            throw new InvalidArgumentException('The statevector must not be empty.');
        }
        $vector = $la->array(array_values($state), dtype: Precision::complexDtype());
        self::assertDimension($vector->shape()[0]);
        return $vector;
    }

    /**
     * |ψ><ψ| を gemm (A A^H) で全成分作る.
     *
     * herk は上三角しか書かないため下三角を PHP ループで補完する必要が
     * あるが、gemm は全成分を書くので補完ループが不要になる。
     * $accumulate を渡すと rho += A A^H として累積する。
     */
    private static function projector(NDArray $column, ?NDArray $accumulate = null): NDArray
    {
        $la = LinearAlgebraProvider::get();
        return $la->gemm(
            $column,
            $column,
            C(1.0, 0.0),
            $accumulate === null ? null : C(1.0, 0.0),
            $accumulate,
            transB: true,
            conjB: true,
        );
    }

    /**
     * |ψ><ψ| を gemm で作る.
     */
    private static function densityFromVector(NDArray $vector): NDArray
    {
        $dim = $vector->shape()[0];
        self::assertDimension($dim);
        if ($dim === 1) {
            $la = LinearAlgebraProvider::get();
            $dtype = $vector->dtype();
            return $la->array([[C(1.0, 0.0)]], dtype: $dtype);
        }
        return self::projector($vector->reshape([$dim, 1]));
    }

    /**
     * 残す qubit を指定して縮約密度行列を直接作る.
     *
     * ψ を [2,2,...,2] に reshape し、残す qubit が行列の行軸に
     * 来るように transpose してから [dA, dB] に reshape し、
     * gemm で A A^H を計算する。全体の |ψ><ψ| を作るより軽い。
     *
     * @param array<int> $keepQubits
     */
    private static function reducedDensityFromVector(NDArray $vector, array $keepQubits): NDArray
    {
        $la = LinearAlgebraProvider::get();
        $dim = $vector->shape()[0];
        self::assertDimension($dim);
        $numQubits = $dim === 1 ? 0 : (int)round(log($dim, 2));
        $keep = self::validateQubits($keepQubits, $numQubits);

        if ($numQubits === 0) {
            return self::densityFromVector($vector);
        }

        $keepSet = array_fill_keys($keep, true);
        $traced = [];
        for ($q = 0; $q < $numQubits; $q++) {
            if (!isset($keepSet[$q])) {
                $traced[] = $q;
            }
        }

        // reshape [2^n] の軸 a は qubit (n-1-a) に対応する。
        // reshape [dA,dB] の行インデックスは上位軸なので、keep の
        // 先頭が最下位ビットになるよう逆順に並べる。
        $rowAxes = [];
        foreach (array_reverse($keep) as $q) {
            $rowAxes[] = $numQubits - 1 - $q;
        }
        $colAxes = [];
        foreach (array_reverse($traced) as $q) {
            $colAxes[] = $numQubits - 1 - $q;
        }
        $perm = array_merge($rowAxes, $colAxes);

        $tensor = $vector->reshape(array_fill(0, $numQubits, 2));
        $permuted = $la->transpose($tensor, $perm);
        $columns = 1 << count($traced);
        $matrix = $permuted->reshape([1 << count($keep), $columns]);
        return self::projector($matrix);
    }

    /**
     * 内部 NDArray を返す. 書き換えないこと.
     */
    public function ndarray(): NDArray
    {
        return $this->ndarray;
    }

    /**
     * @return array<int,array<int,object{real:float,imag:float}>>
     */
    public function data(): array
    {
        return $this->ndarray->toArray();
    }

    /**
     * @return array<int,array<int,object{real:float,imag:float}>>
     */
    public function toArray(): array
    {
        return $this->data();
    }

    public function count(): int
    {
        return $this->dim;
    }

    public function numQubits(): int
    {
        return $this->numQubits;
    }

    /** @deprecated use numQubits() */
    public function num_qubits(): int
    {
        return $this->numQubits();
    }

    public function dim(): int
    {
        return $this->dim;
    }

    public function dimension(): int
    {
        return $this->dim;
    }

    /**
     * 密度行列のトレース Tr(ρ) を返す.
     */
    public function trace(): Complex
    {
        return self::traceOf($this->ndarray);
    }

    /**
     * 純度 Tr(ρ^2) を返す.
     *
     * エルミート行列では Tr(ρ^2) = Σ_ij |ρ_ij|^2 (Frobenius ノルムの二乗)
     * なので、BLAS の nrm2 で計算する。
     */
    public function purity(): float
    {
        $la = LinearAlgebraProvider::get();
        $norm = $la->nrm2($this->ndarray->reshape([$this->dim * $this->dim]));
        return $norm * $norm;
    }

    /**
     * 計算基底の測定確率を返す.
     *
     * 対角成分の実部をそのまま返す (Qiskit と同じ)。数値誤差で
     * 負になった成分だけ 0 に丸める。正規化されていない密度行列では
     * 和が 1 にならないことがある。
     *
     * @param array<int>|null $qargs 部分系のqubitリスト. nullなら全qubit.
     * @return array<int,float> 基底インデックス => 確率
     */
    public function probabilities(?array $qargs = null, ?int $decimals = null): array
    {
        $diagonal = [];
        for ($i = 0; $i < $this->dim; $i++) {
            $value = $this->ndarray[$i][$i];
            $diagonal[$i] = max(0.0, (float)$value->real);
        }

        if ($qargs === null) {
            if ($decimals !== null) {
                foreach ($diagonal as $i => $prob) {
                    $diagonal[$i] = round($prob, $decimals);
                }
            }
            return $diagonal;
        }

        $qargs = self::validateQubits($qargs, $this->numQubits);
        $out = array_fill(0, 1 << count($qargs), 0.0);
        foreach ($diagonal as $index => $prob) {
            $target = 0;
            foreach ($qargs as $position => $q) {
                $target |= ((($index >> $q) & 1) << $position);
            }
            $out[$target] += $prob;
        }
        foreach ($out as $i => $prob) {
            $prob = max(0.0, $prob);
            $out[$i] = $decimals === null ? $prob : round($prob, $decimals);
        }
        return $out;
    }

    /**
     * 計算基底の測定確率をビット列キーの辞書で返す.
     *
     * @param array<int>|null $qargs 部分系のqubitリスト. nullなら全qubit.
     * @return array<int|string,float> ビット列 => 確率 (ゼロ確率は除外)
     */
    public function probabilitiesDict(?array $qargs = null, ?int $decimals = null): array
    {
        $probs = $this->probabilities($qargs, $decimals);
        $bits = $qargs === null ? $this->numQubits : count($qargs);
        $dict = [];
        foreach ($probs as $index => $prob) {
            if ($prob == 0.0) {
                continue;
            }
            $dict[IntMath::int2bit($index, $bits)] = $prob;
        }
        return $dict;
    }

    /** @param array<int>|null $qargs @return array<int|string,float> */
    public function probabilities_dict(?array $qargs = null, ?int $decimals = null): array
    {
        return $this->probabilitiesDict($qargs, $decimals);
    }

    /**
     * <H> = Tr(ρ H) を計算する.
     *
     * SparsePauliOp / Pauli はエルミート行列なので hemm (H * ρ) を使う。
     * 一般の行列は gemm で積を作る。
     *
     * @param SparsePauliOp|Pauli|NDArray|array $oper
     */
    public function expectationValue(SparsePauliOp|Pauli|NDArray|array $oper): float
    {
        $la = LinearAlgebraProvider::get();
        $dtype = $this->ndarray->dtype();
        if ($oper instanceof SparsePauliOp) {
            if ($oper->numQubits() !== $this->numQubits) {
                throw new InvalidArgumentException(
                    'The operator qubit count does not match the density matrix.'
                );
            }
            $product = $la->hemm($oper->toMatrix($dtype), $this->ndarray);
        } elseif ($oper instanceof Pauli) {
            if ($oper->numQubits() !== $this->numQubits) {
                throw new InvalidArgumentException(
                    'The operator qubit count does not match the density matrix.'
                );
            }
            $product = $la->hemm($oper->toMatrix($dtype), $this->ndarray);
        } else {
            if (is_array($oper)) {
                $oper = $la->array($oper, dtype: $dtype);
            } else {
                $operDtype = $oper->dtype();
                if ($operDtype !== NDArray::complex64 && $operDtype !== NDArray::complex128) {
                    $oper = $la->array($oper->toArray(), dtype: $dtype);
                } elseif ($operDtype !== $dtype) {
                    throw new InvalidArgumentException(
                        'The operator precision does not match the density matrix.'
                    );
                }
            }
            $shape = $oper->shape();
            if (count($shape) !== 2 || $shape[0] !== $this->dim || $shape[1] !== $this->dim) {
                throw new InvalidArgumentException(
                    'The operator shape does not match the density matrix.'
                );
            }
            $product = $la->gemm($oper, $this->ndarray);
        }

        $value = self::traceOf($product);
        if (abs($value->imag) > 1.0e-5) {
            throw new InvalidArgumentException(
                "Expectation value has a non-negligible imaginary part: {$value->imag}"
            );
        }
        return (float)$value->real;
    }

    /** @param SparsePauliOp|Pauli|NDArray|array $oper */
    public function expectation_value(SparsePauliOp|Pauli|NDArray|array $oper): float
    {
        return $this->expectationValue($oper);
    }

    /**
     * 指定した qubit をトレースアウトした縮約密度行列を返す.
     *
     * Qiskit の partial_trace(state, qargs) と同じく、$qargs は
     * 「トレースアウトする」qubit のリスト。残る qubit は元の順序
     * (0 から昇順) で並ぶ。
     *
     * @param array<int> $qargs トレースアウトする qubit
     */
    public function partialTrace(array $qargs): self
    {
        $qargs = self::validateQubits($qargs, $this->numQubits);
        if ($qargs === []) {
            return new self($this);
        }
        $traceSet = array_fill_keys($qargs, true);
        $keep = [];
        for ($q = 0; $q < $this->numQubits; $q++) {
            if (!isset($traceSet[$q])) {
                $keep[] = $q;
            }
        }
        return $this->reducedDensityMatrix($keep);
    }

    /** @param array<int> $qargs トレースアウトする qubit */
    public function partial_trace(array $qargs): self
    {
        return $this->partialTrace($qargs);
    }

    /**
     * 残す qubit を指定して縮約密度行列を返す.
     *
     * $qargs の先頭が結果の最下位ビットになる (probabilities() と同じ順序)。
     * 空配列を渡すとトレースを 1x1 密度行列として返す。
     *
     * @param array<int> $qargs 残す qubit
     */
    public function reducedDensityMatrix(array $qargs): self
    {
        $keep = self::validateQubits($qargs, $this->numQubits);
        if ($this->numQubits === 0) {
            return new self($this);
        }
        if (count($keep) === 0) {
            $trace = $this->trace();
            $la = LinearAlgebraProvider::get();
            $matrix = $la->array([[C($trace->real, $trace->imag)]], dtype: $this->ndarray->dtype());
            return new self($matrix);
        }

        $la = LinearAlgebraProvider::get();
        $keepSet = array_fill_keys($keep, true);
        $traced = [];
        for ($q = 0; $q < $this->numQubits; $q++) {
            if (!isset($keepSet[$q])) {
                $traced[] = $q;
            }
        }
        $dimA = 1 << count($keep);
        $dimB = 1 << count($traced);

        // 元の順序をそのまま残す場合はコピーで済ませる。
        if ($keep === range(0, $this->numQubits - 1)) {
            return new self($this);
        }

        // Tr_B(ρ)[a,b] = Σ_t ρ[(a,t),(b,t)] を gemv で計算する。
        //
        // ρ を [2,...,2] (2n 軸) に reshape し、残す軸 (行/列) と
        // トレースする軸 (行/列) の順に transpose して
        // [dA*dA, dB*dB] に reshape する。列方向の単位行列を
        // vec(I_B) に縮約すれば、対角ブロックの和が一度の gemv で得られる。
        $numQubits = $this->numQubits;
        $tensor = $this->ndarray->reshape(array_fill(0, 2 * $numQubits, 2));
        $rowKeep = [];
        $colKeep = [];
        foreach (array_reverse($keep) as $q) {
            $rowKeep[] = $numQubits - 1 - $q;
            $colKeep[] = 2 * $numQubits - 1 - $q;
        }
        $rowTraced = [];
        $colTraced = [];
        foreach (array_reverse($traced) as $q) {
            $rowTraced[] = $numQubits - 1 - $q;
            $colTraced[] = 2 * $numQubits - 1 - $q;
        }
        $perm = array_merge($rowKeep, $colKeep, $rowTraced, $colTraced);
        $permuted = $la->transpose($tensor, $perm);
        $matrix = $permuted->reshape([$dimA * $dimA, $dimB * $dimB]);

        // vec(I_B): 対角位置だけ 1 の選択ベクトル。
        $selector = $la->zeros([$dimB * $dimB], $this->ndarray->dtype());
        for ($t = 0; $t < $dimB; $t++) {
            $selector[$t * $dimB + $t] = C(1.0, 0.0);
        }
        $flat = $la->gemv($matrix, $selector);
        return new self($flat->reshape([$dimA, $dimA]));
    }

    /** @param array<int> $qargs 残す qubit */
    public function reduced_density_matrix(array $qargs): self
    {
        return $this->reducedDensityMatrix($qargs);
    }

    /**
     * エルミート性の簡易チェック.
     *
     * D = H^H - H を BLAS で作り、|D| の最大成分が $atol 以下かを
     * 調べる (Qiskit の is_hermitian と同じ複素絶対値の判定)。
     */
    public function isHermitian(float $atol = 1.0e-8): bool
    {
        $la = LinearAlgebraProvider::get();
        // adj = H^H を作り、その場で H^H - H にする。
        $adj = $la->omatcopy($this->ndarray, trans: true, conj: true);
        $la->axpy($this->ndarray, $adj, C(-1.0, 0.0));
        $flat = $adj->reshape([$this->dim * $this->dim]);
        $index = $la->iamax($flat);
        $value = $flat[$index];
        $modulus = sqrt((float)$value->real ** 2 + (float)$value->imag ** 2);
        return $modulus <= $atol;
    }

    /**
     * 純粋状態の密度行列から状態ベクトルを復元する.
     *
     * 最大対角成分の列を使うので、大域位相は失われる。
     */
    public function toStatevector(): Statevector
    {
        // 純粋状態 ρ = |ψ><ψ| から、最大対角成分の列を使って |ψ> を復元する。
        $purity = $this->purity();
        if (abs($purity - 1.0) > 1.0e-6) {
            throw new InvalidArgumentException('The density matrix is not a pure state.');
        }
        if (!$this->isHermitian()) {
            throw new InvalidArgumentException('The density matrix is not Hermitian.');
        }

        // toArray() を避け、NDArray の要素を直接読む。
        $best = 0;
        $bestValue = -1.0;
        for ($i = 0; $i < $this->dim; $i++) {
            $value = (float)$this->ndarray[$i][$i]->real;
            if ($value > $bestValue) {
                $bestValue = $value;
                $best = $i;
            }
        }
        $scale = sqrt(max($bestValue, 0.0));
        if ($scale === 0.0) {
            throw new InvalidArgumentException('The density matrix has no population.');
        }
        $la = LinearAlgebraProvider::get();
        // 最大対角成分の列を取り出し、sqrt(population) で規格化する。
        // gemv で列を抽出するので PHP ループで Complex を組み立てない。
        $selector = $la->zeros([$this->dim], $this->ndarray->dtype());
        $selector[$best] = C(1.0, 0.0);
        $amplitudes = $la->gemv($this->ndarray, $selector);
        $la->scal(C(1.0 / $scale, 0.0), $amplitudes);
        return new Statevector($amplitudes);
    }

    public function to_statevector(): Statevector
    {
        return $this->toStatevector();
    }

    private static function traceOf(NDArray $matrix): Complex
    {
        $dim = $matrix->shape()[0];
        $real = 0.0;
        $imag = 0.0;
        for ($i = 0; $i < $dim; $i++) {
            $value = $matrix[$i][$i];
            $real += (float)$value->real;
            $imag += (float)$value->imag;
        }
        return new Complex($real, $imag);
    }

    /**
     * @param array<int> $qargs
     * @return array<int>
     */
    private static function validateQubits(array $qargs, int $numQubits): array
    {
        $out = [];
        $seen = [];
        foreach ($qargs as $q) {
            if (!is_int($q) || $q < 0 || $q >= $numQubits) {
                throw new InvalidArgumentException("Qubit index out of range: {$q}");
            }
            if (isset($seen[$q])) {
                throw new InvalidArgumentException("Duplicated qargs: {$q}");
            }
            $seen[$q] = true;
            $out[] = $q;
        }
        return $out;
    }
}
