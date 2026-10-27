<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Math;

use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Math\NelderMead;

class NelderMeadTest extends TestCase
{
    public function testMinimizeQuadratic(): void
    {
        [$x, $f] = NelderMead::minimize(
            static fn(array $p): float => ($p[0] - 1.0) ** 2 + ($p[1] + 2.0) ** 2,
            [0.0, 0.0],
            maxIter: 500,
            tol: 1.0e-10
        );

        $this->assertEqualsWithDelta(1.0, $x[0], 1.0e-4);
        $this->assertEqualsWithDelta(-2.0, $x[1], 1.0e-4);
        $this->assertLessThan(1.0e-8, $f);
    }

    public function testMinimizeOneDimensionalFunction(): void
    {
        [$x, $f] = NelderMead::minimize(
            static fn(array $p): float => ($p[0] - 3.0) ** 2,
            [0.0],
            step: 0.5,
            maxIter: 300,
            tol: 1.0e-10
        );

        $this->assertEqualsWithDelta(3.0, $x[0], 1.0e-4);
        $this->assertLessThan(1.0e-8, $f);
    }

    public function testMaxIterationsIsRespected(): void
    {
        $calls = 0;
        NelderMead::minimize(
            static function (array $p) use (&$calls): float {
                $calls++;
                return ($p[0] - 100.0) ** 2;
            },
            [0.0],
            maxIter: 3
        );

        $this->assertGreaterThan(0, $calls);
    }
}
