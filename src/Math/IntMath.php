<?php
declare(strict_types=1);

namespace Rindow\Qubit\Math;

class IntMath
{
    public static function gcd(int $a, int $b): int
    {
        $a = abs($a);
        $b = abs($b);
        while ($b != 0) {
            [$a, $b] = [$b, $a % $b];
        }
        return $a;
    }

    public static function powmod(int $base, int $exp, int $mod): int
    {
        $r = 1 % $mod;
        $base %= $mod;
        while ($exp > 0) {
            if ($exp & 1) {
                $r = ($r * $base) % $mod;
            }
            $base = ($base * $base) % $mod;
            $exp >>= 1;
        }
        return $r;
    }
}