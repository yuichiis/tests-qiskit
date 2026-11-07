<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Math;

use InvalidArgumentException;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Math\Minimize;
use Rindow\Qubit\Math\OptimizeResult;

class MinimizeTest extends TestCase
{
    public function testCobylaReturnsOptimizeResult(): void
    {
        $result = Minimize::minimize(
            static fn(array $p): float => ($p[0] - 2.0) ** 2 + ($p[1] - 3.0) ** 2,
            [0.0, 0.0],
            method: 'COBYLA',
            options: ['maxiter' => 500, 'tol' => 1.0e-8],
        );

        $this->assertInstanceOf(OptimizeResult::class, $result);
        $this->assertIsArray($result->x);
        $this->assertEqualsWithDelta(2.0, $result->x[0], 1.0e-3);
        $this->assertEqualsWithDelta(3.0, $result->x[1], 1.0e-3);
        $this->assertEqualsWithDelta(0.0, $result->fun, 1.0e-5);
        $this->assertTrue($result->success);
        $this->assertSame(0, $result->status);
        $this->assertSame('Optimization terminated successfully.', $result->message);
        $this->assertGreaterThan(0, $result->nfev);
    }

    public function testCobylaIsDefaultMethod(): void
    {
        $result = Minimize::minimize(
            static fn(array $p): float => ($p[0] + 1.0) ** 2,
            [0.0],
        );

        $this->assertEqualsWithDelta(-1.0, $result->x[0], 1.0e-3);
    }

    public function testMethodNameIsCaseInsensitive(): void
    {
        $result = Minimize::minimize(
            static fn(array $p): float => ($p[0] - 4.0) ** 2,
            [0.0],
            method: 'cobyla',
            options: ['maxiter' => 300],
        );

        $this->assertEqualsWithDelta(4.0, $result->x[0], 1.0e-3);
    }

    public function testConstraintsAreForwarded(): void
    {
        $result = Minimize::minimize(
            static fn(array $p): float => $p[0] ** 2,
            [-2.0],
            method: 'COBYLA',
            options: ['maxiter' => 500, 'tol' => 1.0e-8],
            constraints: [static fn(array $p): float => $p[0] - 1.0],
        );

        $this->assertGreaterThanOrEqual(1.0 - 1.0e-3, $result->x[0]);
    }

    public function testNelderMeadReturnsOptimizeResult(): void
    {
        $result = Minimize::minimize(
            static fn(array $p): float => ($p[0] - 5.0) ** 2,
            [0.0],
            method: 'Nelder-Mead',
            options: ['maxiter' => 300, 'tol' => 1.0e-10, 'step' => 0.5],
        );

        $this->assertInstanceOf(OptimizeResult::class, $result);
        $this->assertEqualsWithDelta(5.0, $result->x[0], 1.0e-4);
        $this->assertEqualsWithDelta(0.0, $result->fun, 1.0e-8);
        $this->assertTrue($result->success);
    }

    public function testNelderMeadUnderscoreVariant(): void
    {
        $result = Minimize::minimize(
            static fn(array $p): float => ($p[0] + 2.0) ** 2,
            [0.0],
            method: 'nelder_mead',
            options: ['maxiter' => 300],
        );

        $this->assertEqualsWithDelta(-2.0, $result->x[0], 1.0e-3);
    }

    public function testUnsupportedMethod(): void
    {
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Unsupported optimization method: BFGS');
        Minimize::minimize(static fn(array $p): float => 0.0, [0.0], method: 'BFGS');
    }

    public function testOptimizeResultIsMutable(): void
    {
        $result = new OptimizeResult([1.0, 2.0], 0.5, 10, true, 0, 'done');

        $this->assertSame([1.0, 2.0], $result->x);
        $this->assertSame(0.5, $result->fun);
        $this->assertSame(10, $result->nfev);
        $this->assertTrue($result->success);
        $this->assertSame(0, $result->status);
        $this->assertSame('done', $result->message);

        $result->fun = 0.1;
        $this->assertSame(0.1, $result->fun);
    }
}
