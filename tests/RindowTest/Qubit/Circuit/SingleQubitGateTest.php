<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Circuit;

use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Circuit\HGate;
use Rindow\Qubit\Circuit\IdGate;
use Rindow\Qubit\Circuit\IGate;
use Rindow\Qubit\Circuit\PhaseGate;
use Rindow\Qubit\Circuit\RXGate;
use Rindow\Qubit\Circuit\RYGate;
use Rindow\Qubit\Circuit\RZGate;
use Rindow\Qubit\Circuit\SdgGate;
use Rindow\Qubit\Circuit\SGate;
use Rindow\Qubit\Circuit\TdgGate;
use Rindow\Qubit\Circuit\TGate;
use Rindow\Qubit\Circuit\UGate;
use Rindow\Qubit\Circuit\XGate;
use Rindow\Qubit\Circuit\YGate;
use Rindow\Qubit\Circuit\ZGate;
use RindowTest\Qubit\TestUtils;

class SingleQubitGateTest extends TestCase
{
    use TestUtils;

    /** @param array<mixed> $args @param array<array<mixed>> $expected */
    #[DataProvider('matrixProvider')]
    public function testMatrix(string $class, array $args, string $name, array $expected): void
    {
        $gate = new $class(...$args);

        $this->assertSame($name, $gate->name());
        $this->assertSame(2, $gate->dimension());
        $this->assertGateMatrix($gate, $expected);
    }

    /** @return array<string,array{class-string,array<mixed>,string,array<array<mixed>>}> */
    public static function matrixProvider(): array
    {
        $s = 1.0 / sqrt(2.0);
        $theta = M_PI / 2.0;
        $phi = M_PI / 3.0;
        $lambda = M_PI / 4.0;
        $c = cos($theta / 2.0);
        $sn = sin($theta / 2.0);

        return [
            'H' => [HGate::class, [], 'h', [[$s, $s], [$s, -$s]]],
            'X' => [XGate::class, [], 'x', [[0, 1], [1, 0]]],
            'Y' => [
                YGate::class,
                [],
                'y',
                [[0, self::c(0.0, -1.0)], [self::c(0.0, 1.0), 0]],
            ],
            'Z' => [ZGate::class, [], 'z', [[1, 0], [0, -1]]],
            'S' => [SGate::class, [], 's', [[1, 0], [0, self::c(0.0, 1.0)]]],
            'Sdg' => [SdgGate::class, [], 'sdg', [[1, 0], [0, self::c(0.0, -1.0)]]],
            'T' => [TGate::class, [], 't', [[1, 0], [0, self::c($s, $s)]]],
            'Tdg' => [TdgGate::class, [], 'tdg', [[1, 0], [0, self::c($s, -$s)]]],
            'I' => [IGate::class, [], 'id', [[1, 0], [0, 1]]],
            'Id' => [IdGate::class, [], 'id', [[1, 0], [0, 1]]],
            'RX' => [
                RXGate::class,
                [$theta],
                'rx',
                [
                    [self::c($c, 0.0), self::c(0.0, -$sn)],
                    [self::c(0.0, -$sn), self::c($c, 0.0)],
                ],
            ],
            'RY' => [
                RYGate::class,
                [$theta],
                'ry',
                [
                    [self::c($c, 0.0), self::c(-$sn, 0.0)],
                    [self::c($sn, 0.0), self::c($c, 0.0)],
                ],
            ],
            'RZ' => [
                RZGate::class,
                [$theta],
                'rz',
                [
                    [self::c(cos($theta / 2.0), -sin($theta / 2.0)), self::c(0.0, 0.0)],
                    [self::c(0.0, 0.0), self::c(cos($theta / 2.0), sin($theta / 2.0))],
                ],
            ],
            'P' => [
                PhaseGate::class,
                [$lambda],
                'p',
                [[1, 0], [0, self::c(cos($lambda), sin($lambda))]],
            ],
            'U' => [
                UGate::class,
                [$theta, $phi, $lambda],
                'u',
                [
                    [self::c($c, 0.0), self::c(-$sn * cos($lambda), -$sn * sin($lambda))],
                    [
                        self::c($sn * cos($phi), $sn * sin($phi)),
                        self::c($c * cos($phi + $lambda), $c * sin($phi + $lambda)),
                    ],
                ],
            ],
        ];
    }

    public function testParameterAccessors(): void
    {
        $rx = new RXGate(0.5);
        $this->assertSame(0.5, $rx->theta());

        $ry = new RYGate(0.25);
        $this->assertSame(0.25, $ry->theta());

        $rz = new RZGate(0.75);
        $this->assertSame(0.75, $rz->theta());

        $p = new PhaseGate(1.25);
        $this->assertSame(1.25, $p->lambda());

        $u = new UGate(0.1, 0.2, 0.3);
        $this->assertSame(0.1, $u->theta());
        $this->assertSame(0.2, $u->phi());
        $this->assertSame(0.3, $u->lambda());
    }

    public function testInverses(): void
    {
        $this->assertInstanceOf(HGate::class, (new HGate())->inverse());
        $this->assertInstanceOf(XGate::class, (new XGate())->inverse());
        $this->assertInstanceOf(YGate::class, (new YGate())->inverse());
        $this->assertInstanceOf(ZGate::class, (new ZGate())->inverse());
        $this->assertInstanceOf(IGate::class, (new IGate())->inverse());

        $this->assertInstanceOf(SdgGate::class, (new SGate())->inverse());
        $this->assertInstanceOf(SGate::class, (new SdgGate())->inverse());
        $this->assertInstanceOf(TdgGate::class, (new TGate())->inverse());
        $this->assertInstanceOf(TGate::class, (new TdgGate())->inverse());

        $this->assertEqualsWithDelta(-0.5, (new RXGate(0.5))->inverse()->theta(), 1.0e-12);
        $this->assertEqualsWithDelta(-0.5, (new RYGate(0.5))->inverse()->theta(), 1.0e-12);
        $this->assertEqualsWithDelta(-0.5, (new RZGate(0.5))->inverse()->theta(), 1.0e-12);
        $this->assertEqualsWithDelta(-0.5, (new PhaseGate(0.5))->inverse()->lambda(), 1.0e-12);

        $u = new UGate(0.1, 0.2, 0.3);
        $inverse = $u->inverse();
        $this->assertInstanceOf(UGate::class, $inverse);
        $this->assertEqualsWithDelta(-0.1, $inverse->theta(), 1.0e-12);
        $this->assertEqualsWithDelta(-0.3, $inverse->phi(), 1.0e-12);
        $this->assertEqualsWithDelta(-0.2, $inverse->lambda(), 1.0e-12);
    }

    public function testDraw(): void
    {
        $gate = new HGate();
        $this->assertSame(
            ['cells' => [1 => 'h'], 'span' => [1]],
            $gate->draw([1])
        );

        $gate = new RXGate(0.5);
        $this->assertSame(
            ['cells' => [0 => 'rx(0.5)'], 'span' => [0]],
            $gate->draw([0])
        );

        $gate = new UGate(0.1, 0.2, 0.3);
        $this->assertSame(
            ['cells' => [2 => 'u(0.1,0.2,0.3)'], 'span' => [2]],
            $gate->draw([2])
        );
    }

    /** @return array{float,float} */
    private static function c(float $real, float $imag): array
    {
        return ['real' => $real, 'imag' => $imag];
    }
}
