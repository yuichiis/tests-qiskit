<?php
declare(strict_types=1);

namespace Rindow\Qubit\Primitives;

/**
 * qiskitの BitArray に相当. `getCounts()` / `get_counts()` でカウント辞書を返す.
 */
class MeasData
{
    /** @var array<string,int> */
    private array $counts;
    private int $numShots;

    /** @param array<string,int> $counts */
    public function __construct(array $counts, int $numShots)
    {
        $this->counts = $counts;
        $this->numShots = $numShots;
    }

    /** @return array<string,int> ビット列 => カウント (降順ソートせず、キー昇順) */
    public function getCounts(): array
    {
        return $this->counts;
    }

    /** qiskit互換エイリアス */
    public function get_counts(): array
    {
        return $this->counts;
    }

    public function numShots(): int
    {
        return $this->numShots;
    }

    /** @return array<string,float> ビット列 => 確率 */
    public function getProbabilities(): array
    {
        $probs = [];
        foreach ($this->counts as $bits => $c) {
            $probs[$bits] = $c / $this->numShots;
        }
        return $probs;
    }

    /** qiskit互換エイリアス */
    public function get_probabilities(): array
    {
        return $this->getProbabilities();
    }
}
