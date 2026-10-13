<?php
declare(strict_types=1);

namespace Rindow\Qubit\Visualization;

/**
 * qiskit.visualization.plot_histogram のPHP版.
 * ターミナル向けASCIIヒストグラムを返す (+echoするヘルパー付き).
 */
class Visualization
{
    /**
     * @param array<string,int> $counts ビット列 => カウント
     */
    public static function plotHistogram(array $counts, int $width = 40): string
    {
        if (empty($counts)) {
            return "(empty counts)\n";
        }
        ksort($counts);
        $max = max($counts);
        $total = array_sum($counts);
        $lines = [];
        foreach ($counts as $bits => $c) {
            $barLen = $max > 0 ? (int)round($c / $max * $width) : 0;
            $bar = str_repeat('█', $barLen);
            $ratio = $total > 0 ? $c / $total : 0.0;
            $lines[] = sprintf('%s : %5d (%5.1f%%) %s', $bits, $c, $ratio * 100.0, $bar);
        }
        return implode("\n", $lines) . "\n";
    }

    /** qiskit互換エイリアス (snake_case) */
    public static function plot_histogram(array $counts, int $width = 40): string
    {
        return self::plotHistogram($counts, $width);
    }

    /** @param array<string,int> $counts */
    public static function showHistogram(array $counts, int $width = 40): void
    {
        echo self::plotHistogram($counts, $width);
    }

    /** qiskit互換エイリアス */
    public static function show_histogram(array $counts, int $width = 40): void
    {
        self::showHistogram($counts, $width);
    }
}
