<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Circuit;

use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Circuit\CPhaseGate;
use Rindow\Qubit\Circuit\CXGate;
use Rindow\Qubit\Circuit\CZGate;
use Rindow\Qubit\Circuit\SwapGate;
use RindowTest\Qubit\TestUtils;

class TwoQubitGateTest extends TestCase
{
    use TestUtils;

    /** @param array<array<mixed>> $expected */
    #[DataProvider('matrixProvider')]
    public function testMatrix(string $class, float $lambda, string $name, array $expected): void
    {
        $gate = $class === CPhaseGate::class
            ? new CPhaseGate($lambda)
            : new $class();

        $this->assertSame($name, $gate->name());
        $this->assertSame(4, $gate->dimension());
        $this->assertGateMatrix($gate, $expected);
    }

    /** @return array<string,array{class-string,float,string,array<array<mixed>>}> */
    public static function matrixProvider(): array
    {
        $lambda = M_PI / 3.0;

        return [
            'CX' => [
                CXGate::class,
                0.0,
                'cx',
                [
                    [1, 0, 0, 0],
                    [0, 0, 0, 1],
                    [0, 0, 1, 0],
                    [0, 1, 0, 0],
                ],
            ],
            'CZ' => [
                CZGate::class,
                0.0,
                'cz',
                [
                    [1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, -1],
                ],
            ],
            'CP' => [
                CPhaseGate::class,
                $lambda,
                'cp',
                [
                    [1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, ['real' => cos($lambda), 'imag' => sin($lambda)]],
                ],
            ],
            'Swap' => [
                SwapGate::class,
                0.0,
                'swap',
                [
                    [1, 0, 0, 0],
                    [0, 0, 1, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 1],
                ],
            ],
        ];
    }

    public function testCPhaseParameterAndInverse(): void
    {
        $gate = new CPhaseGate(0.5);
        $this->assertSame(0.5, $gate->lambda());

        $inverse = $gate->inverse();
        $this->assertInstanceOf(CPhaseGate::class, $inverse);
        $this->assertEqualsWithDelta(-0.5, $inverse->lambda(), 1.0e-12);
    }

    public function testInvolutoryInverses(): void
    {
        $this->assertInstanceOf(CXGate::class, (new CXGate())->inverse());
        $this->assertInstanceOf(CZGate::class, (new CZGate())->inverse());

        $swap = new SwapGate();
        $inverse = $swap->inverse();
        $this->assertInstanceOf(SwapGate::class, $inverse);
        $this->assertGateMatrix($inverse, [
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 1, 0, 0],
            [0, 0, 0, 1],
        ]);
    }

    public function testDraw(): void
    {
        $this->assertSame(
            ['cells' => [0 => '●', 1 => '⊕'], 'span' => [0, 1]],
            (new CXGate())->draw([0, 1])
        );
        $this->assertSame(
            ['cells' => [0 => '●', 1 => '●'], 'span' => [0, 1]],
            (new CZGate())->draw([0, 1])
        );
        $this->assertSame(
            ['cells' => [0 => '●', 1 => 'cp(0.5)'], 'span' => [0, 1]],
            (new CPhaseGate(0.5))->draw([0, 1])
        );
        $this->assertSame(
            ['cells' => [0 => '×', 1 => '×'], 'span' => [0, 1]],
            (new SwapGate())->draw([0, 1])
        );
    }

    public function testDrawWithIntermediateQubit(): void
    {
        $this->assertEquals(
            ['cells' => [0 => '●', 2 => '⊕', 1 => '│'], 'span' => [0, 1, 2]],
            (new CXGate())->draw([0, 2])
        );
    }
}
