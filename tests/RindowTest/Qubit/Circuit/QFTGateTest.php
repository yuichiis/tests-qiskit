<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Circuit;

use InvalidArgumentException;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Circuit\QFTGate;
use RindowTest\Qubit\TestUtils;

class QFTGateTest extends TestCase
{
    use TestUtils;

    public function testOneQubitMatrix(): void
    {
        $gate = new QFTGate(1);
        $s = 1.0 / sqrt(2.0);

        $this->assertSame('qft', $gate->name());
        $this->assertSame(1, $gate->numQubits());
        $this->assertSame(2, $gate->dimension());
        $this->assertFalse($gate->isInverse());
        $this->assertGateMatrix($gate, [[$s, $s], [$s, -$s]]);
    }

    public function testTwoQubitMatrix(): void
    {
        $gate = new QFTGate(2);
        $half = 0.5;

        $this->assertSame(2, $gate->numQubits());
        $this->assertSame(4, $gate->dimension());
        $this->assertGateMatrix($gate, [
            [self::c($half), self::c($half), self::c($half), self::c($half)],
            [self::c($half), self::c(0.0, $half), self::c(-$half), self::c(0.0, -$half)],
            [self::c($half), self::c(-$half), self::c($half), self::c(-$half)],
            [self::c($half), self::c(0.0, -$half), self::c(-$half), self::c(0.0, $half)],
        ]);
    }

    public function testInverse(): void
    {
        $forward = new QFTGate(2);
        $inverse = $forward->inverse();

        $this->assertInstanceOf(QFTGate::class, $inverse);
        $this->assertTrue($inverse->isInverse());
        $this->assertSame('qft_dg', $inverse->name());
        $this->assertSame(2, $inverse->numQubits());

        $forwardValues = $forward->parameters();
        $inverseValues = $inverse->parameters();
        foreach ($forwardValues as $i => $value) {
            $this->assertEqualsWithDelta($value->real, $inverseValues[$i]->real, 1.0e-6);
            $this->assertEqualsWithDelta(-$value->imag, $inverseValues[$i]->imag, 1.0e-6);
        }
    }

    public function testDraw(): void
    {
        $gate = new QFTGate(2);
        $this->assertSame(
            ['cells' => [0 => 'QFT', 1 => 'QFT'], 'span' => [0, 1]],
            $gate->draw([0, 1])
        );
        $this->assertEquals(
            ['cells' => [0 => 'QFT', 2 => 'QFT', 1 => '│'], 'span' => [0, 1, 2]],
            $gate->draw([0, 2])
        );
        $this->assertSame(
            ['cells' => [0 => 'QFT†', 1 => 'QFT†'], 'span' => [0, 1]],
            $gate->inverse()->draw([0, 1])
        );
    }

    public function testInvalidNumQubits(): void
    {
        $this->expectException(InvalidArgumentException::class);
        new QFTGate(0);
    }

    public function testTooManyQubits(): void
    {
        $this->expectException(InvalidArgumentException::class);
        new QFTGate(31);
    }

    /** @return array{float,float} */
    private static function c(float $real, float $imag = 0.0): array
    {
        return ['real' => $real, 'imag' => $imag];
    }
}
