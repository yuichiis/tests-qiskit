<?php
declare(strict_types=1);

namespace Rindow\Qubit;

use InvalidArgumentException;

trait CircuitUtils
{
    // ---- 単一量子ビットゲート ----
    // Qiskit互換: qubit引数は int または intの配列/Traversableを受け付け、
    // 配列の場合は各qubitに順にゲートを適用する.

    /**
     * @param int|array<int>|\Traversable $qubit
     * @return array<int> 正規化されたqubitリスト
     */
    private static function normalizeQubitList(int|array|\Traversable $qubit): array
    {
        if ($qubit instanceof \Traversable) {
            $qubit = iterator_to_array($qubit, false);
        }
        if (is_array($qubit)) {
            $list = array_values($qubit);
            if (count($list) === 0) {
                throw new InvalidArgumentException('Qubit list must not be empty.');
            }
            return $list;
        }
        return [$qubit];
    }

    /**
     * @param int|array<int>|\Traversable $clbit
     * @return array<int>
     */
    private static function normalizeClbitList(int|array|\Traversable $clbit): array
    {
        if ($clbit instanceof \Traversable) {
            $clbit = iterator_to_array($clbit, false);
        }
        if (is_array($clbit)) {
            $list = array_values($clbit);
            if (count($list) === 0) {
                throw new InvalidArgumentException('Clbit list must not be empty.');
            }
            return $list;
        }
        return [$clbit];
    }

    /**
     * 2qubitゲート用の引数展開. Qiskit互換のブロードキャスト:
     * - 両方単一 => 1ペア
     * - 片方が単一、もう一方がリスト => 単一を各要素にブロードキャスト
     * - 両方リストで同数 => 1対1対応
     * - それ以外 => エラー
     * @return array<array{int,int}>
     */
    private static function expandTwoQubitArgs(int|array|\Traversable $q1, int|array|\Traversable $q2): array
    {
        $list1 = self::normalizeQubitList($q1);
        $list2 = self::normalizeQubitList($q2);
        $n1 = count($list1);
        $n2 = count($list2);
        if ($n1 === 1 && $n2 === 1) {
            return [[(int)$list1[0], (int)$list2[0]]];
        }
        if ($n1 === 1) {
            $pairs = [];
            foreach ($list2 as $b) {
                $pairs[] = [(int)$list1[0], (int)$b];
            }
            return $pairs;
        }
        if ($n2 === 1) {
            $pairs = [];
            foreach ($list1 as $a) {
                $pairs[] = [(int)$a, (int)$list2[0]];
            }
            return $pairs;
        }
        if ($n1 !== $n2) {
            throw new InvalidArgumentException(
                "Qubit list size mismatch: {$n1} vs {$n2}."
            );
        }
        $pairs = [];
        foreach ($list1 as $i => $a) {
            $pairs[] = [(int)$a, (int)$list2[$i]];
        }
        return $pairs;
    }
}