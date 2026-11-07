<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Math;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Circuit\HGate;
use Rindow\Qubit\Math\Precision;
use Rindow\Qubit\Quantum\StatevectorSimulator;
use Rindow\Qubit\QuantumInfo\Pauli;
use Rindow\Qubit\QuantumInfo\SparsePauliOp;
use Rindow\Qubit\QuantumInfo\Statevector;
use function Rindow\Math\Matrix\C;

class PrecisionTest extends TestCase
{
    protected function tearDown(): void
    {
        Precision::reset();
    }

    public function testDefaultIsFloat64AndComplex128(): void
    {
        Precision::reset();

        $this->assertSame(NDArray::float64, Precision::realDtype());
        $this->assertSame(NDArray::complex128, Precision::complexDtype());
    }

    public function testSetFromRealDtype(): void
    {
        Precision::set(NDArray::float32);
        $this->assertSame(NDArray::float32, Precision::realDtype());
        $this->assertSame(NDArray::complex64, Precision::complexDtype());

        Precision::set(NDArray::float64);
        $this->assertSame(NDArray::float64, Precision::realDtype());
        $this->assertSame(NDArray::complex128, Precision::complexDtype());
    }

    public function testSetFromComplexDtype(): void
    {
        Precision::set(NDArray::complex64);
        $this->assertSame(NDArray::float32, Precision::realDtype());
        $this->assertSame(NDArray::complex64, Precision::complexDtype());

        Precision::set(NDArray::complex128);
        $this->assertSame(NDArray::float64, Precision::realDtype());
        $this->assertSame(NDArray::complex128, Precision::complexDtype());
    }

    public function testSetMatchingPair(): void
    {
        Precision::set(NDArray::float32, NDArray::complex64);
        $this->assertSame(NDArray::complex64, Precision::complexDtype());

        Precision::set(NDArray::float64, NDArray::complex128);
        $this->assertSame(NDArray::complex128, Precision::complexDtype());
    }

    public function testMismatchedPairIsRejected(): void
    {
        $this->expectException(InvalidArgumentException::class);
        Precision::set(NDArray::float32, NDArray::complex128);
    }

    public function testInvalidDtypeIsRejected(): void
    {
        $this->expectException(InvalidArgumentException::class);
        Precision::set(NDArray::int32);
    }

    public function testAllComponentsFollowPrecision(): void
    {
        foreach ([
            [NDArray::float32, NDArray::complex64],
            [NDArray::float64, NDArray::complex128],
        ] as [$realDtype, $complexDtype]) {
            Precision::set($realDtype, $complexDtype);

            $this->assertSame($complexDtype, (new HGate())->matrix()->dtype());
            $this->assertSame(
                $complexDtype,
                (new StatevectorSimulator(1))->stateNDArray()->dtype()
            );
            $this->assertSame($complexDtype, (new Pauli('X'))->toMatrix()->dtype());
            $this->assertSame(
                $complexDtype,
                SparsePauliOp::fromList([['Z', 1.0]])->toMatrix()->dtype()
            );
            $this->assertSame(
                $complexDtype,
                (new Statevector([C(1.0, 0.0), C(0.0, 0.0)]))->ndarray()->dtype()
            );
        }
    }
}
