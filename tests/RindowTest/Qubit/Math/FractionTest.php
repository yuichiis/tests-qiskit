<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Math;

use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Math\Fraction;

class FractionTest extends TestCase
{
    /** @param array{int,int} $expected */
    #[DataProvider('fractionProvider')]
    public function testLimitDenominator(int $num, int $den, int $maxDen, array $expected): void
    {
        $this->assertSame($expected, Fraction::limitDenominator($num, $den, $maxDen));
    }

    /** @return array<string,array{int,int,int,array{int,int}}> */
    public static function fractionProvider(): array
    {
        return [
            'pi' => [355, 113, 100, [22, 7]],
            'half' => [1, 2, 10, [1, 2]],
            'two quarters' => [2, 4, 10, [1, 2]],
            'integer' => [7, 1, 10, [7, 1]],
            'zero' => [0, 1, 10, [0, 1]],
        ];
    }

    public function testLimitDenominatorRespectsMaximumDenominator(): void
    {
        [$p, $q] = Fraction::limitDenominator(355, 113, 10);
        $this->assertLessThanOrEqual(10, $q);
        $this->assertSame(22, $p);
        $this->assertSame(7, $q);
    }
}
