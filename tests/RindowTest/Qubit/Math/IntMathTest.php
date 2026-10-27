<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Math;

use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Math\IntMath;

class IntMathTest extends TestCase
{
    public function testGcd(): void
    {
        $this->assertSame(6, IntMath::gcd(12, 18));
        $this->assertSame(6, IntMath::gcd(18, 12));
        $this->assertSame(6, IntMath::gcd(-12, 18));
        $this->assertSame(6, IntMath::gcd(12, -18));
        $this->assertSame(5, IntMath::gcd(0, 5));
        $this->assertSame(5, IntMath::gcd(5, 0));
        $this->assertSame(0, IntMath::gcd(0, 0));
        $this->assertSame(1, IntMath::gcd(17, 13));
    }

    public function testInt2Bit(): void
    {
        $this->assertSame('0101', IntMath::int2bit(5, 4));
        $this->assertSame('1010', IntMath::int2bit(10, 4));
        $this->assertSame('0', IntMath::int2bit(0, 1));
        $this->assertSame('00000000', IntMath::int2bit(0, 8));
        $this->assertSame('1', IntMath::int2bit(1, 0));
        $this->assertSame('101', IntMath::int2bit(5, 3));
    }

    public function testPowmod(): void
    {
        $this->assertSame(24, IntMath::powmod(2, 10, 1000));
        $this->assertSame(1, IntMath::powmod(2, 0, 7));
        $this->assertSame(8, IntMath::powmod(5, 3, 13));
        $this->assertSame(4, IntMath::powmod(10, 2, 12));
        $this->assertSame(0, IntMath::powmod(6, 3, 9));
    }
}
