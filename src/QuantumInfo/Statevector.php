<?php
declare(strict_types=1);

namespace Rindow\Qubit\QuantumInfo;

use Countable;
use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use IteratorAggregate;
use Rindow\Qubit\Math\IntMath;
use Rindow\Qubit\Math\LinearAlgebraProvider;
use Rindow\Qubit\Math\Precision;
use Rindow\Qubit\Quantum\StatevectorSimulator;
use Rindow\Qubit\QuantumCircuit;
use Traversable;
use function Rindow\Math\Matrix\C;

/**
 * qiskit.quantum_info.Statevector のPHP版.
 *
 * 内部では NDArray を保持し、シミュレータから受け取った状態ベクトルを
 * PHP の float に戻さずに BLAS 演算へ渡せる。PHP 配列が必要なときだけ
 * data() / getStatevector() で遅延変換する。
 *
 * ```php
 * use Rindow\Qubit\QuantumInfo\Statevector;
 * $probs = (new Statevector($qc))->probabilitiesDict();
 * ```
 */
class Statevector implements IteratorAggregate, Countable
{
    /** @var array<int,object{real:float,imag:float}>|null 遅延変換キャッシュ */
    private ?array $statevector = null;
    private ?NDArray $ndarray = null;
    private int $numQubits;
    private int $dim;

    /**
     * @param QuantumCircuit|StatevectorSimulator|Statevector|NDArray|array<int,object{real:float,imag:float}> $data
     */
    public function __construct(
        QuantumCircuit|StatevectorSimulator|Statevector|NDArray|array $data
    ) {
        if ($data instanceof self) {
            $this->statevector = $data->statevector;
            $this->ndarray = $data->ndarray;
            $this->numQubits = $data->numQubits;
            $this->dim = $data->dim;
            return;
        }
        if ($data instanceof QuantumCircuit) {
            $sim = new StatevectorSimulator($data->numQubits());
            $sim->applyCircuit($data);
            $this->ndarray = $this->copyNDArray($sim->stateNDArray());
            $this->numQubits = $data->numQubits();
            $this->dim = 1 << $this->numQubits;
            return;
        }
        if ($data instanceof StatevectorSimulator) {
            $this->ndarray = $this->copyNDArray($data->stateNDArray());
            $this->numQubits = $data->numQubits();
            $this->dim = $data->dimension();
            return;
        }
        if ($data instanceof NDArray) {
            $shape = $data->shape();
            if (count($shape) !== 1) {
                throw new InvalidArgumentException(
                    'The statevector NDArray must be one-dimensional.'
                );
            }
            $dim = $shape[0];
            $this->assertDimension($dim);
            $this->ndarray = $this->copyNDArray($data);
            $this->dim = $dim;
            $this->numQubits = (int)round(log($dim, 2));
            return;
        }
        $dim = count($data);
        $this->assertDimension($dim);
        $this->statevector = array_values($data);
        $this->dim = $dim;
        $this->numQubits = (int)round(log($dim, 2));
    }

    private function assertDimension(int $dim): void
    {
        if ($dim < 1 || ($dim & ($dim - 1)) !== 0) {
            throw new InvalidArgumentException(
                'The length of the statevector must be a power of two.'
            );
        }
    }

    private function copyNDArray(NDArray $state): NDArray
    {
        $la = LinearAlgebraProvider::get();
        if ($state->dtype() !== NDArray::complex64
            && $state->dtype() !== NDArray::complex128) {
            throw new InvalidArgumentException(
                'The statevector NDArray must have a complex dtype.'
            );
        }
        return $la->copy($state);
    }

    /**
     * 内部 NDArray を返す. 変換が必要な場合だけ PHP 配列から構築する.
     */
    public function ndarray(): NDArray
    {
        if ($this->ndarray !== null) {
            return $this->ndarray;
        }
        $la = LinearAlgebraProvider::get();
        $values = [];
        foreach ($this->statevector as $amp) {
            $values[] = C(
                (float)($amp->real ?? 0.0),
                (float)($amp->imag ?? 0.0),
            );
        }
        $this->ndarray = $la->array($values, dtype: Precision::complexDtype());
        return $this->ndarray;
    }

    /**
     * @return array<int,object{real:float,imag:float}> 各要素は real/imag を持つ複素数オブジェクト
     */
    public function getStatevector(): array
    {
        return $this->toPhpArray();
    }

    /**
     * @return array<int,object{real:float,imag:float}>
     */
    public function data(): array
    {
        return $this->toPhpArray();
    }

    /**
     * @return array<int,object{real:float,imag:float}>
     */
    private function toPhpArray(): array
    {
        if ($this->statevector !== null) {
            return $this->statevector;
        }
        $la = LinearAlgebraProvider::get();
        $flat = $this->ndarray->reshape([$this->dim]);
        $out = [];
        for ($i = 0; $i < $this->dim; $i++) {
            $value = $la->scalar($flat[$i]);
            $out[] = (object)[
                'real' => (float)$value->real,
                'imag' => (float)$value->imag,
            ];
        }
        $this->statevector = $out;
        return $out;
    }

    /**
     * @return Traversable<int,object{real:float,imag:float}>
     */
    public function getIterator(): Traversable
    {
        return new \ArrayIterator($this->toPhpArray());
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
     * 計算基底の測定確率を返す.
     *
     * @param array<int>|null $qargs 部分系のqubitリスト. nullなら全qubit.
     * @return array<int,float> 基底インデックス => 確率
     */
    public function probabilities(?array $qargs = null, ?int $decimals = null): array
    {
        $la = LinearAlgebraProvider::get();
        // |ψ|^2 は NDArray 内で計算し、PHP へは確率配列だけを取り出す。
        $probs = $la->square($la->abs($this->ndarray()))
            ->reshape([$this->dim])
            ->toArray();

        if ($qargs === null) {
            foreach ($probs as $i => $prob) {
                $prob = min(1.0, max(0.0, (float)$prob));
                $probs[$i] = $decimals === null ? $prob : round($prob, $decimals);
            }
            return $probs;
        }

        $qargs = array_values($qargs);
        if (count($qargs) > $this->numQubits) {
            throw new InvalidArgumentException('Too many qargs for the statevector.');
        }
        $seen = [];
        foreach ($qargs as $q) {
            if (!is_int($q) || $q < 0 || $q >= $this->numQubits) {
                throw new InvalidArgumentException("Qubit index out of range: {$q}");
            }
            if (isset($seen[$q])) {
                throw new InvalidArgumentException("Duplicated qargs: {$q}");
            }
            $seen[$q] = true;
        }

        $out = array_fill(0, 1 << count($qargs), 0.0);
        foreach ($probs as $index => $prob) {
            $target = 0;
            foreach ($qargs as $position => $q) {
                $target |= ((($index >> $q) & 1) << $position);
            }
            $out[$target] += (float)$prob;
        }
        foreach ($out as $i => $prob) {
            $prob = min(1.0, max(0.0, $prob));
            $out[$i] = $decimals === null ? $prob : round($prob, $decimals);
        }
        return $out;
    }

    /**
     * 計算基底の測定確率をビット列キーの辞書で返す.
     *
     * PHPの配列キーは先頭に0がない数値文字列 (例: "1000") がintに
     * 変換されるため、必要に応じて利用側で (string) キャストする.
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
}
