<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Math;

use PHPUnit\Framework\TestCase;
use Rindow\Math\Matrix\MatrixOperator;
use Rindow\Qubit\Math\LinearAlgebraProvider;

class LinearAlgebraProviderTest extends TestCase
{
    public function testSharedInstance(): void
    {
        $first = LinearAlgebraProvider::get();
        $second = LinearAlgebraProvider::get();

        $this->assertSame($first, $second);
    }

    public function testSetAndRestore(): void
    {
        $original = LinearAlgebraProvider::get();
        $mo = new MatrixOperator();
        try {
            LinearAlgebraProvider::set($mo->la());
            $this->assertSame($mo->la(), LinearAlgebraProvider::get());
        } finally {
            LinearAlgebraProvider::set($original);
        }

        $this->assertSame($original, LinearAlgebraProvider::get());
    }

    public function testGetLinearAlgebraWithoutAcceleration(): void
    {
        $env = getenv('RINDOW_QUBIT_LINEAR_ALGEBRA');
        putenv('RINDOW_QUBIT_LINEAR_ALGEBRA');
        try {
            $mo = new MatrixOperator();
            $this->assertSame($mo->la(), LinearAlgebraProvider::getLinearAlgebra($mo));
        } finally {
            if ($env === false) {
                putenv('RINDOW_QUBIT_LINEAR_ALGEBRA');
            } else {
                putenv('RINDOW_QUBIT_LINEAR_ALGEBRA=' . $env);
            }
        }
    }

    public function testInvalidAccelerationService(): void
    {
        $this->expectException(\RuntimeException::class);
        LinearAlgebraProvider::getLinearAlgebra(new MatrixOperator(), 'unknown');
    }
}
