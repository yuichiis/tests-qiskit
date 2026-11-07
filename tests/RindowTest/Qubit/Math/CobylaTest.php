<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Math;

use InvalidArgumentException;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Math\Cobyla;

class CobylaTest extends TestCase
{
    public function testMinimizeQuadratic(): void
    {
        $nfev = 0;
        $fun = static function (array $x) use (&$nfev): float {
            $nfev++;
            return ($x[0] - 1.0) ** 2 + ($x[1] + 2.0) ** 2;
        };

        [$x, $f, $evaluations] = Cobyla::minimize($fun, [0.0, 0.0], maxiter: 500, tol: 1.0e-8);

        $this->assertEqualsWithDelta(1.0, $x[0], 1.0e-3);
        $this->assertEqualsWithDelta(-2.0, $x[1], 1.0e-3);
        $this->assertEqualsWithDelta(0.0, $f, 1.0e-5);
        $this->assertSame($nfev, $evaluations);
        $this->assertGreaterThan(0, $evaluations);
    }

    public function testOneDimensionalQuadratic(): void
    {
        [$x, $f] = Cobyla::minimize(
            static fn(array $p): float => ($p[0] - 3.0) ** 2,
            [0.0],
            maxiter: 300,
            tol: 1.0e-10,
        );

        $this->assertEqualsWithDelta(3.0, $x[0], 1.0e-5);
        $this->assertEqualsWithDelta(0.0, $f, 1.0e-8);
    }

    public function testMaxIterationsLimitsEvaluations(): void
    {
        $fun = static fn(array $p): float => ($p[0] - 10.0) ** 2;
        [$x, $f, $nfev] = Cobyla::minimize($fun, [0.0], maxiter: 1, tol: 1.0e-12);

        // 初回評価 + 勾配プローブ + 候補評価
        $this->assertSame(3, $nfev);
        $this->assertLessThan(10.0, $x[0]);
        $this->assertGreaterThan(0.0, $f);
    }

    public function testZeroMaxIterations(): void
    {
        [$x, $f, $nfev] = Cobyla::minimize(
            static fn(array $p): float => ($p[0] - 2.0) ** 2,
            [5.0],
            maxiter: 0,
        );

        $this->assertSame([5.0], $x);
        $this->assertSame(9.0, $f);
        $this->assertSame(1, $nfev);
    }

    public function testConstraintIsSatisfied(): void
    {
        // x >= 1 の制約下で (x - 0)^2 を最小化 → x = 1 が期待値
        $constraint = static fn(array $p): float => $p[0] - 1.0;
        [$x, $f] = Cobyla::minimize(
            static fn(array $p): float => $p[0] ** 2,
            [-3.0],
            [$constraint],
            maxiter: 500,
            tol: 1.0e-8,
        );

        $this->assertGreaterThanOrEqual(1.0 - 1.0e-3, $x[0]);
        $this->assertEqualsWithDelta(1.0, $x[0], 5.0e-3);
        $this->assertEqualsWithDelta(1.0, $f, 1.0e-3);
    }

    public function testEmptyInitialPoint(): void
    {
        $this->expectException(InvalidArgumentException::class);
        Cobyla::minimize(static fn(array $p): float => 0.0, []);
    }

    public function testNegativeMaxIterations(): void
    {
        $this->expectException(InvalidArgumentException::class);
        Cobyla::minimize(static fn(array $p): float => 0.0, [0.0], maxiter: -1);
    }

    public function testInvalidTolerance(): void
    {
        $this->expectException(InvalidArgumentException::class);
        Cobyla::minimize(static fn(array $p): float => 0.0, [0.0], tol: 0.0);
    }

    public function testInvalidRhobeg(): void
    {
        $this->expectException(InvalidArgumentException::class);
        Cobyla::minimize(static fn(array $p): float => 0.0, [0.0], rhobeg: -1.0);
    }
}
