<?php
declare(strict_types=1);

namespace Rindow\Qubit\QuantumInfo;

use Countable;
use InvalidArgumentException;
use IteratorAggregate;
use Rindow\Qubit\Math\IntMath;
use Rindow\Qubit\Quantum\StatevectorSimulator;
use Rindow\Qubit\QuantumCircuit;
use Traversable;

/**
 * qiskit.quantum_info.Statevector のPHP版.
 *
 * ```php
 * use Rindow\Qubit\QuantumInfo\Statevector;
 * $probs = (new Statevector($qc))->probabilitiesDict();
 * ```
 */
class Statevector implements IteratorAggregate, Countable
{
    /** @var array<int,object{real:float,imag:float}> */
    private array $statevector;
    private int $numQubits;
    private int $dim;

    /**
     * @param QuantumCircuit|StatevectorSimulator|Statevector|array<int,object{real:float,imag:float}> $data
     */
    public function __construct(QuantumCircuit|StatevectorSimulator|Statevector|array $data)
    {
        if ($data instanceof self) {
            $this->statevector = $data->statevector;
            $this->numQubits = $data->numQubits;
            $this->dim = $data->dim;
            return;
        }
        if ($data instanceof QuantumCircuit) {
            $sim = new StatevectorSimulator($data->numQubits());
            $sim->applyCircuit($data);
            $this->statevector = $sim->getStatevector();
            $this->numQubits = $data->numQubits();
            $this->dim = 1 << $this->numQubits;
            return;
        }
        if ($data instanceof StatevectorSimulator) {
            $this->statevector = $data->getStatevector();
            $this->numQubits = $data->numQubits();
            $this->dim = $data->dimension();
            return;
        }
        $dim = count($data);
        if ($dim < 1 || ($dim & ($dim - 1)) !== 0) {
            throw new InvalidArgumentException(
                'The length of the statevector must be a power of two.'
            );
        }
        $this->statevector = array_values($data);
        $this->dim = $dim;
        $this->numQubits = (int)round(log($dim, 2));
    }

    /**
     * @return array<int,object{real:float,imag:float}> 各要素は real/imag を持つ複素数オブジェクト
     */
    public function getStatevector(): array
    {
        return $this->statevector;
    }

    /**
     * @return array<int,object{real:float,imag:float}>
     */
    public function data(): array
    {
        return $this->statevector;
    }

    /**
     * @return Traversable<int,object{real:float,imag:float}>
     */
    public function getIterator(): Traversable
    {
        return new \ArrayIterator($this->statevector);
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
        if ($qargs === null) {
            $probs = [];
            foreach ($this->statevector as $amp) {
                $prob = $this->amplitudeProbability($amp);
                $prob = min(1.0, max(0.0, $prob));
                $probs[] = $decimals === null ? $prob : round($prob, $decimals);
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

        $probs = array_fill(0, 1 << count($qargs), 0.0);
        foreach ($this->statevector as $index => $amp) {
            $target = 0;
            foreach ($qargs as $position => $q) {
                $target |= ((($index >> $q) & 1) << $position);
            }
            $probs[$target] += $this->amplitudeProbability($amp);
        }
        foreach ($probs as $i => $prob) {
            $prob = min(1.0, max(0.0, $prob));
            $probs[$i] = $decimals === null ? $prob : round($prob, $decimals);
        }
        return $probs;
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

    /**
     * @param object{real:float,imag:float} $amp
     */
    private function amplitudeProbability(object $amp): float
    {
        $real = (float)$amp->real;
        $imag = (float)$amp->imag;
        return $real * $real + $imag * $imag;
    }
}
