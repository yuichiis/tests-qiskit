<?php
declare(strict_types=1);

namespace Rindow\Qubit\Math;

class Fraction
{
    /** 連分数で分母上限つきの近似分数 [p, q] を返す */
    public static function limitDenominator(int $num, int $den, int $maxDen): array
    {
        $cf = [];
        while ($den != 0) {
            $cf[] = intdiv($num, $den);
            [$num, $den] = [$den, $num % $den];
        }
        $p0 = 1;
        $p1 = $cf[0];
        $q0 = 0;
        $q1 = 1;
        $best = [$p1, $q1];
        for ($i = 1; $i < count($cf); $i++) {
            $p2 = $cf[$i] * $p1 + $p0;
            $q2 = $cf[$i] * $q1 + $q0;
            if ($q2 > $maxDen) {
                break;
            }
            $best = [$p2, $q2];
            [$p0, $p1] = [$p1, $p2];
            [$q0, $q1] = [$q1, $q2];
        }
        return $best;
    }
}
