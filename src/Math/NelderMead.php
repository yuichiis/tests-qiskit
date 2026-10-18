<?php
declare(strict_types=1);

namespace Rindow\Qubit\Math;

class NelderMead
{
    /**
     * Nelder-Mead法（微分不要の古典オプティマイザ）
     * minimize function f in N dimensions.
     * 
     * @param callable $f 目的関数
     * @param array $x0 初期値
     * @param float $step 初期ステップ幅
     * @param int $maxIter 最大反復回数
     * @param float $tol 収束許容誤差
     * @return array [$best_x, $best_f]
     */
    public static function minimize(callable $f, array $x0, float $step = 0.3, int $maxIter = 200, float $tol = 1e-6): array
    {
        $n = count($x0);
        $simplex = [$x0];
        for ($i = 0; $i < $n; $i++) {
            $q = $x0;
            $q[$i] += $step;
            $simplex[] = $q;
        }
        $fvals = array_map($f, $simplex);
        for ($iter = 0; $iter < $maxIter; $iter++) {
            array_multisort($fvals, SORT_ASC, SORT_NUMERIC, $simplex);
            $maxd = 0.0;
            foreach ($simplex as $q) {
                $d = 0.0;
                foreach ($q as $j => $v) {
                    $d += ($v - $simplex[0][$j]) ** 2;
                }
                $maxd = max($maxd, $d);
            }
            if (sqrt($maxd) < $tol) {
                break;
            }
            $cent = array_fill(0, $n, 0.0);
            for ($i = 0; $i < $n; $i++) {
                foreach ($simplex[$i] as $j => $v) {
                    $cent[$j] += $v / $n;
                }
            }
            $mirror = static function (array $c, array $w, float $a) {
                $out = [];
                foreach ($c as $j => $v) {
                    $out[$j] = $v + $a * ($v - $w[$j]);
                }
                return $out;
            };
            $worst = $simplex[$n];
            $refl = $mirror($cent, $worst, 1.0);
            $fr = $f($refl);
            if ($fr < $fvals[0]) {
                $exp = $mirror($cent, $worst, 2.0);
                $fe = $f($exp);
                if ($fe < $fr) {
                    $simplex[$n] = $exp;
                    $fvals[$n] = $fe;
                } else {
                    $simplex[$n] = $refl;
                    $fvals[$n] = $fr;
                }
            } elseif ($fr < $fvals[$n - 1]) {
                $simplex[$n] = $refl;
                $fvals[$n] = $fr;
            } else {
                [$tmp, $ft] = $fr < $fvals[$n] ? [$refl, $fr] : [$worst, $fvals[$n]];
                $con = $mirror($cent, $tmp, 0.5);
                $fc = $f($con);
                if ($fc < $ft) {
                    $simplex[$n] = $con;
                    $fvals[$n] = $fc;
                } else {
                    for ($i = 1; $i <= $n; $i++) {
                        foreach ($simplex[$i] as $j => $v) {
                            $simplex[$i][$j] = $simplex[0][$j] + 0.5 * ($v - $simplex[0][$j]);
                        }
                        $fvals[$i] = $f($simplex[$i]);
                    }
                }
            }
        }
        array_multisort($fvals, SORT_ASC, SORT_NUMERIC, $simplex);
        return [$simplex[0], $fvals[0]];
    }
}