<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Math;

use InvalidArgumentException;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Math\HermitianEigen;
use Rindow\Qubit\QuantumInfo\SparsePauliOp;

class HermitianEigenTest extends TestCase
{
    public function testOneByOne(): void
    {
        $this->assertSame([3.5], HermitianEigen::eigvalsh([[3.5]]));
    }

    public function testTwoByTwoDiagonal(): void
    {
        $eigenvalues = HermitianEigen::eigvalsh([[2, 0], [0, -1]]);

        $this->assertEqualsWithDelta([-1.0, 2.0], $eigenvalues, 1.0e-12);
    }

    public function testTwoByTwoSymmetric(): void
    {
        // [[2, 1], [1, 2]] の固有値は 1, 3
        $eigenvalues = HermitianEigen::eigvalsh([[2, 1], [1, 2]]);

        $this->assertEqualsWithDelta([1.0, 3.0], $eigenvalues, 1.0e-12);
    }

    public function testThreeByThree(): void
    {
        // 2x2 ブロック [[2,-1],[-1,2]] の固有値 1, 3 と、対角成分 2
        $matrix = [
            [2.0, -1.0, 0.0],
            [-1.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
        ];
        $eigenvalues = HermitianEigen::eigvalsh($matrix);

        $this->assertEqualsWithDelta([1.0, 2.0, 3.0], $eigenvalues, 1.0e-12);
    }

    public function testH2HamiltonianMatchesQiskit(): void
    {
        $h = SparsePauliOp::fromList([
            ['II', -1.052373245772859],
            ['IZ', 0.39793742484318045],
            ['ZI', -0.39793742484318045],
            ['ZZ', -0.01128010425623538],
            ['XX', 0.18093119978423156],
        ]);

        $eigenvalues = HermitianEigen::eigvalsh(
            $h->toMatrix(\Interop\Polite\Math\Matrix\NDArray::complex128)
        );

        $this->assertEqualsWithDelta([
            -1.8572750302023795,
            -1.2445845498133257,
            -0.8827221502448627,
            -0.2249112528308669,
        ], $eigenvalues, 1.0e-10);
    }

    public function testH2HamiltonianComplex64Matrix(): void
    {
        $h = SparsePauliOp::fromList([
            ['II', -1.052373245772859],
            ['IZ', 0.39793742484318045],
            ['ZI', -0.39793742484318045],
            ['ZZ', -0.01128010425623538],
            ['XX', 0.18093119978423156],
        ]);

        $matrix = $h->toMatrix(\Interop\Polite\Math\Matrix\NDArray::complex64);

        $this->assertSame(\Interop\Polite\Math\Matrix\NDArray::complex64, $matrix->dtype());
        $eigenvalues = HermitianEigen::eigvalsh($matrix);
        $this->assertEqualsWithDelta(-1.8572750302023795, $eigenvalues[0], 1.0e-6);
    }

    public function testSortedAscending(): void
    {
        $eigenvalues = HermitianEigen::eigvalsh([
            [4.0, 1.0, 0.0],
            [1.0, 4.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);

        $sorted = $eigenvalues;
        sort($sorted, SORT_NUMERIC);
        $this->assertSame($sorted, $eigenvalues);
    }

    public function testSlightAsymmetryIsAveraged(): void
    {
        $eigenvalues = HermitianEigen::eigvalsh([[2.0, 1.0], [1.0 + 1.0e-12, 2.0]]);

        $this->assertEqualsWithDelta([1.0, 3.0], $eigenvalues, 1.0e-10);
    }

    public function testAcceptsNDArray(): void
    {
        $la = \Rindow\Qubit\Math\LinearAlgebraProvider::get();
        $matrix = $la->array([[2.0, 1.0], [1.0, 2.0]], dtype: \Interop\Polite\Math\Matrix\NDArray::float64);

        $this->assertEqualsWithDelta([1.0, 3.0], HermitianEigen::eigvalsh($matrix), 1.0e-12);
    }

    public function testNonSquareMatrix(): void
    {
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix must be square.');
        HermitianEigen::eigvalsh([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
    }

    public function testNonSymmetricMatrix(): void
    {
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix must be symmetric.');
        HermitianEigen::eigvalsh([[1.0, 2.0], [3.0, 4.0]]);
    }

    public function testComplexMatrixRejected(): void
    {
        $this->expectException(InvalidArgumentException::class);
        HermitianEigen::eigvalsh([
            [[1.0, 0.0], [0.0, 0.5]],
            [[0.0, -0.5], [1.0, 0.0]],
        ]);
    }

    public function testEmptyMatrix(): void
    {
        $this->expectException(InvalidArgumentException::class);
        HermitianEigen::eigvalsh([]);
    }
}
