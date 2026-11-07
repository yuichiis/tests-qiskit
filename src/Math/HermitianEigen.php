<?php
declare(strict_types=1);

namespace Rindow\Qubit\Math;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Math\Matrix\Complex;

/**
 * エルミート行列の固有値を求める小さなソルバ.
 *
 * numpy.linalg.eigvalsh の代替として、実対称行列 (および実数値の
 * エルミート行列) にヤコビ回転法を適用する。VQE の厳密対角化のような
 * 小さな行列を想定している。
 *
 * ```php
 * $eigenvalues = HermitianEigen::eigvalsh($hamiltonian->toMatrix());
 * $groundEnergy = min($eigenvalues);
 * ```
 */
final class HermitianEigen
{
    /**
     * @param array<int,array<int,mixed>>|NDArray $matrix
     * @return array<int,float> 昇順の固有値
     */
    public static function eigvalsh(array|NDArray $matrix): array
    {
        $a = self::toRealMatrix($matrix);
        $n = count($a);
        if ($n === 0) {
            throw new InvalidArgumentException('Matrix must not be empty.');
        }
        if ($n === 1) {
            return [$a[0][0]];
        }

        $scale = 0.0;
        for ($i = 0; $i < $n; $i++) {
            $scale = max($scale, abs($a[$i][$i]));
        }
        $eps = 1.0e-14 * max(1.0, $scale);

        // ヤコビ回転法。非対角成分が十分小さくなるまで掃き出す。
        for ($sweep = 0; $sweep < 100; $sweep++) {
            $off = 0.0;
            for ($p = 0; $p < $n - 1; $p++) {
                for ($q = $p + 1; $q < $n; $q++) {
                    $off += $a[$p][$q] * $a[$p][$q];
                }
            }
            if (sqrt($off) <= $eps) {
                break;
            }
            for ($p = 0; $p < $n - 1; $p++) {
                for ($q = $p + 1; $q < $n; $q++) {
                    $apq = $a[$p][$q];
                    if (abs($apq) <= 1.0e-18) {
                        continue;
                    }
                    $app = $a[$p][$p];
                    $aqq = $a[$q][$q];
                    $tau = ($aqq - $app) / (2.0 * $apq);
                    if ($tau >= 0.0) {
                        $t = 1.0 / ($tau + sqrt(1.0 + $tau * $tau));
                    } else {
                        $t = -1.0 / (-$tau + sqrt(1.0 + $tau * $tau));
                    }
                    $c = 1.0 / sqrt(1.0 + $t * $t);
                    $s = $t * $c;

                    $a[$p][$p] = $app - $t * $apq;
                    $a[$q][$q] = $aqq + $t * $apq;
                    $a[$p][$q] = 0.0;
                    $a[$q][$p] = 0.0;

                    for ($k = 0; $k < $n; $k++) {
                        if ($k === $p || $k === $q) {
                            continue;
                        }
                        $akp = $a[$k][$p];
                        $akq = $a[$k][$q];
                        $newKp = $c * $akp - $s * $akq;
                        $newKq = $s * $akp + $c * $akq;
                        $a[$k][$p] = $newKp;
                        $a[$p][$k] = $newKp;
                        $a[$k][$q] = $newKq;
                        $a[$q][$k] = $newKq;
                    }
                }
            }
        }

        $eigenvalues = [];
        for ($i = 0; $i < $n; $i++) {
            $eigenvalues[] = $a[$i][$i];
        }
        sort($eigenvalues, SORT_NUMERIC);
        return $eigenvalues;
    }

    /**
     * @param array<int,array<int,mixed>>|NDArray $matrix
     * @return array<int,array<int,float>>
     */
    private static function toRealMatrix(array|NDArray $matrix): array
    {
        if ($matrix instanceof NDArray) {
            $shape = $matrix->shape();
            if (count($shape) !== 2 || $shape[0] !== $shape[1]) {
                throw new InvalidArgumentException('Matrix must be square.');
            }
            $matrix = $matrix->toArray();
        }
        $n = count($matrix);
        if ($n === 0) {
            throw new InvalidArgumentException('Matrix must not be empty.');
        }
        $out = [];
        foreach ($matrix as $i => $row) {
            if (!is_array($row) || count($row) !== $n) {
                throw new InvalidArgumentException('Matrix must be square.');
            }
            $values = array_values($row);
            $outRow = [];
            foreach ($values as $j => $value) {
                [$real, $imag] = self::complexParts($value);
                if (abs($imag) > 1.0e-9) {
                    throw new InvalidArgumentException(
                        'Only real symmetric matrices are supported for eigenvalue computation.'
                    );
                }
                $outRow[] = $real;
            }
            $out[] = $outRow;
        }
        for ($i = 0; $i < $n; $i++) {
            for ($j = $i + 1; $j < $n; $j++) {
                $diff = abs($out[$i][$j] - $out[$j][$i]);
                if ($diff > 1.0e-9 * max(1.0, abs($out[$i][$j]), abs($out[$j][$i]))) {
                    throw new InvalidArgumentException('Matrix must be symmetric.');
                }
                $average = 0.5 * ($out[$i][$j] + $out[$j][$i]);
                $out[$i][$j] = $average;
                $out[$j][$i] = $average;
            }
        }
        return $out;
    }

    /** @return array{float,float} */
    private static function complexParts(mixed $value): array
    {
        if (is_int($value) || is_float($value)) {
            return [(float)$value, 0.0];
        }
        if ($value instanceof Complex) {
            return [(float)$value->real, (float)$value->imag];
        }
        if (is_object($value) && isset($value->real)) {
            return [(float)$value->real, (float)($value->imag ?? 0.0)];
        }
        if (is_array($value)) {
            return [
                (float)($value['real'] ?? $value[0] ?? 0.0),
                (float)($value['imag'] ?? $value[1] ?? 0.0),
            ];
        }
        throw new InvalidArgumentException('Unsupported matrix element type.');
    }
}
