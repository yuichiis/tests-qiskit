<?php
declare(strict_types=1);

namespace RindowTest\Qubit\QuantumInfo;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use PHPUnit\Framework\TestCase;
use Rindow\Math\Matrix\Complex;
use Rindow\Qubit\Math\LinearAlgebraProvider;
use Rindow\Qubit\Math\Precision;
use Rindow\Qubit\QuantumInfo\DensityMatrix;
use Rindow\Qubit\QuantumInfo\SparsePauliOp;
use Rindow\Qubit\QuantumInfo\Statevector;
use Rindow\Qubit\QuantumCircuit;
use function Rindow\Math\Matrix\C;

class DensityMatrixTest extends TestCase
{
    private function h2(): SparsePauliOp
    {
        return SparsePauliOp::fromList([
            ['II', -1.052373245772859],
            ['IZ', 0.39793742484318045],
            ['ZI', -0.39793742484318045],
            ['ZZ', -0.01128010425623538],
            ['XX', 0.18093119978423156],
        ]);
    }

    private function bellStatevector(): Statevector
    {
        $qc = (new QuantumCircuit(2))->h(0)->cx(0, 1);
        return new Statevector($qc);
    }

    public function testFromComputationalBasis(): void
    {
        $rho = new DensityMatrix([C(1.0, 0.0), C(0.0, 0.0)]);

        $this->assertSame(1, $rho->numQubits());
        $this->assertSame(2, $rho->dim());
        $this->assertSame(2, $rho->dimension());
        $this->assertCount(2, $rho);
        $this->assertEqualsWithDelta(1.0, $rho->trace()->real, 1.0e-12);
        $this->assertEqualsWithDelta(0.0, $rho->trace()->imag, 1.0e-12);
        $this->assertEqualsWithDelta(1.0, $rho->purity(), 1.0e-12);
        $this->assertFloatVectorEquals([1.0, 0.0], $rho->probabilities());
        $this->assertComplexMatrix([[1, 0], [0, 0]], $rho->data());
        $this->assertSame(['0' => 1.0], $rho->probabilitiesDict());
    }

    public function testFromCircuitBellState(): void
    {
        $rho = new DensityMatrix(
            (new QuantumCircuit(2))->h(0)->cx(0, 1)
        );

        $this->assertSame(2, $rho->numQubits());
        $this->assertSame(4, $rho->dim());
        $this->assertComplexMatrix([
            [0.5, 0.0, 0.0, 0.5],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0, 0.5],
        ], $rho->data());
        $this->assertEqualsWithDelta(1.0, $rho->trace()->real, 1.0e-12);
        $this->assertEqualsWithDelta(1.0, $rho->purity(), 1.0e-12);
        $this->assertFloatVectorEquals([0.5, 0.0, 0.0, 0.5], $rho->probabilities());
        $this->assertTrue($rho->isHermitian());
    }

    public function testFromStatevectorMatchesOuterProduct(): void
    {
        $amplitudes = [
            C(0.5, 0.1),
            C(-0.2, 0.3),
            C(0.4, -0.25),
            C(0.0, 0.6),
        ];
        $la = LinearAlgebraProvider::get();
        $vector = $la->array($amplitudes, dtype: NDArray::complex128);
        $rho = DensityMatrix::fromStatevector($vector);

        $this->assertSame(NDArray::complex128, $rho->ndarray()->dtype());
        $this->assertSame([4, 4], $rho->ndarray()->shape());

        $expected = [];
        for ($i = 0; $i < 4; $i++) {
            $row = [];
            for ($j = 0; $j < 4; $j++) {
                $a = $amplitudes[$i];
                $b = $amplitudes[$j];
                $row[] = C(
                    $a->real * $b->real + $a->imag * $b->imag,
                    $a->imag * $b->real - $a->real * $b->imag,
                );
            }
            $expected[] = $row;
        }
        $this->assertComplexMatrix($expected, $rho->data(), 1.0e-12);
        // 入力ベクトルは規格化されていないので、トレースはノルムの二乗。
        $this->assertEqualsWithDelta(0.9725, $rho->trace()->real, 1.0e-12);
        $this->assertEqualsWithDelta(0.9725 * 0.9725, $rho->purity(), 1.0e-12);
    }

    public function testFromStatevectorWithComplex64(): void
    {
        $la = LinearAlgebraProvider::get();
        $vector = $la->array([
            C(0.5, 0.0),
            C(0.0, 0.5),
            C(0.5, 0.0),
            C(0.0, -0.5),
        ], dtype: NDArray::complex64);

        $rho = DensityMatrix::fromStatevector($vector);

        $this->assertSame(NDArray::complex64, $rho->ndarray()->dtype());
        $this->assertEqualsWithDelta(1.0, $rho->trace()->real, 1.0e-5);
        $this->assertEqualsWithDelta(1.0, $rho->purity(), 1.0e-5);
    }

    public function testFromEnsembleNormalizes(): void
    {
        $rho = DensityMatrix::fromEnsemble([
            [C(1.0, 0.0), C(0.0, 0.0)],
            [C(0.0, 0.0), C(1.0, 0.0)],
        ]);

        $this->assertComplexMatrix([
            [0.5, 0.0],
            [0.0, 0.5],
        ], $rho->data(), 1.0e-12);
        $this->assertEqualsWithDelta(1.0, $rho->trace()->real, 1.0e-12);
        $this->assertEqualsWithDelta(0.5, $rho->purity(), 1.0e-12);
        $this->assertFloatVectorEquals([0.5, 0.5], $rho->probabilities());
    }

    public function testFromEnsembleWithoutNormalization(): void
    {
        $rho = DensityMatrix::fromEnsemble([
            [C(1.0, 0.0), C(0.0, 0.0)],
            [C(0.0, 0.0), C(1.0, 0.0)],
        ], normalize: false);

        $this->assertComplexMatrix([
            [1.0, 0.0],
            [0.0, 1.0],
        ], $rho->data(), 1.0e-12);
        $this->assertEqualsWithDelta(2.0, $rho->trace()->real, 1.0e-12);
        $this->assertEqualsWithDelta(2.0, $rho->purity(), 1.0e-12);
        // 正規化していないので確率の和はトレース (2) になる。
        $this->assertFloatVectorEquals([1.0, 1.0], $rho->probabilities());
    }

    public function testFromEnsembleMixedPrecisionThrows(): void
    {
        $la = LinearAlgebraProvider::get();
        $a = $la->array([C(1.0, 0.0), C(0.0, 0.0)], dtype: NDArray::complex64);
        $b = $la->array([C(0.0, 0.0), C(0.0, 1.0)], dtype: NDArray::complex128);

        $this->expectException(InvalidArgumentException::class);
        DensityMatrix::fromEnsemble([$a, $b]);
    }

    public function testFromEnsembleMixedPrecisionThrowsInAnyOrder(): void
    {
        $la = LinearAlgebraProvider::get();
        $a = $la->array([C(1.0, 0.0), C(0.0, 0.0)], dtype: NDArray::complex64);
        $b = $la->array([C(0.0, 0.0), C(0.0, 1.0)], dtype: NDArray::complex128);

        $this->expectException(InvalidArgumentException::class);
        DensityMatrix::fromEnsemble([$b, $a]);
    }

    public function testFromInt(): void
    {
        $rho = DensityMatrix::fromInt(1, 2);

        $this->assertSame(2, $rho->numQubits());
        $this->assertComplexMatrix([
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ], $rho->data());
        $this->assertFloatVectorEquals([0.0, 1.0, 0.0, 0.0], $rho->probabilities());
    }

    public function testPartialTraceOfBellState(): void
    {
        $rho = new DensityMatrix($this->bellStatevector());

        $this->assertComplexMatrix([
            [0.5, 0.0],
            [0.0, 0.5],
        ], $rho->partialTrace([1])->data(), 1.0e-12);
        $this->assertComplexMatrix([
            [0.5, 0.0],
            [0.0, 0.5],
        ], $rho->partialTrace([0])->data(), 1.0e-12);

        // 何もトレースアウトしない場合はコピー。
        $copy = $rho->partialTrace([]);
        $this->assertComplexMatrix($rho->data(), $copy->data(), 1.0e-12);

        // すべてトレースアウトすると 1x1 のトレースになる。
        $scalar = $rho->partialTrace([0, 1]);
        $this->assertSame(1, $scalar->dim());
        $this->assertSame(0, $scalar->numQubits());
        $this->assertEqualsWithDelta(1.0, $scalar->data()[0][0]->real, 1.0e-12);
    }

    public function testReducedDensityMatrixOrdering(): void
    {
        // |+>_q0 |0>_q1. 状態ベクトルの基底は q0 が LSB。
        $rho = new DensityMatrix((new QuantumCircuit(2))->h(0));

        $this->assertComplexMatrix([
            [0.5, 0.5],
            [0.5, 0.5],
        ], $rho->reducedDensityMatrix([0])->data(), 1.0e-12);
        $this->assertComplexMatrix([
            [1.0, 0.0],
            [0.0, 0.0],
        ], $rho->reducedDensityMatrix([1])->data(), 1.0e-12);

        // qargs の先頭が結果の LSB になる。
        $this->assertComplexMatrix([
            [0.5, 0.0, 0.5, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.5, 0.0, 0.5, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ], $rho->reducedDensityMatrix([1, 0])->data(), 1.0e-12);

        // 元の順序を保つ全 qubit 指定は元の行列と同じ。
        $this->assertComplexMatrix(
            $rho->data(),
            $rho->reducedDensityMatrix([0, 1])->data(),
            1.0e-12
        );
    }

    public function testExpectationValueMatchesSparsePauliOp(): void
    {
        $h = $this->h2();
        $state = $this->bellStatevector();
        $rho = new DensityMatrix($state);

        $expected = $h->expectationValue($state);
        $this->assertEqualsWithDelta($expected, $rho->expectationValue($h), 1.0e-12);
        $this->assertEqualsWithDelta($expected, $h->expectationValue($rho), 1.0e-12);
        $this->assertEqualsWithDelta(
            $expected,
            $rho->expectationValue($h->toMatrix()),
            1.0e-12
        );
        $this->assertEqualsWithDelta(
            $expected,
            $rho->expectation_value($h),
            1.0e-12
        );
    }

    public function testExpectationValueOnBasisState(): void
    {
        $rho = DensityMatrix::fromInt(0, 2);

        $this->assertEqualsWithDelta(
            -1.0636533500290943,
            $rho->expectationValue($this->h2()),
            1.0e-6
        );
    }

    public function testToStatevectorFromPureState(): void
    {
        $rho = new DensityMatrix($this->bellStatevector());
        $state = $rho->toStatevector();

        $this->assertSame(2, $state->numQubits());
        $this->assertFloatVectorEquals([0.5, 0.0, 0.0, 0.5], $state->probabilities());
        $this->assertEquals($state->data(), $rho->to_statevector()->data());
    }

    public function testToStatevectorRejectsMixedState(): void
    {
        $rho = DensityMatrix::fromEnsemble([
            [C(1.0, 0.0), C(0.0, 0.0)],
            [C(0.0, 0.0), C(1.0, 0.0)],
        ]);

        $this->expectException(InvalidArgumentException::class);
        $rho->toStatevector();
    }

    public function testIsHermitianDetectsNonHermitianMatrix(): void
    {
        $la = LinearAlgebraProvider::get();

        $nonConjugate = new DensityMatrix($la->array([
            [C(1.0, 0.0), C(0.0, 0.0)],
            [C(1.0, 0.0), C(0.0, 0.0)],
        ], dtype: NDArray::complex128));
        $this->assertFalse($nonConjugate->isHermitian());

        $imaginaryDiagonal = new DensityMatrix($la->array([
            [C(1.0, 0.1), C(0.0, 0.0)],
            [C(0.0, 0.0), C(0.0, 0.0)],
        ], dtype: NDArray::complex128));
        $this->assertFalse($imaginaryDiagonal->isHermitian());
    }

    public function testExpectationValueRejectsMixedPrecisionOperator(): void
    {
        $la = LinearAlgebraProvider::get();
        $rho = new DensityMatrix($la->array([
            [C(1.0, 0.0), C(0.0, 0.0)],
            [C(0.0, 0.0), C(0.0, 0.0)],
        ], dtype: NDArray::complex128));
        $oper = $la->array([
            [C(1.0, 0.0), C(0.0, 0.0)],
            [C(0.0, 0.0), C(1.0, 0.0)],
        ], dtype: NDArray::complex64);

        $this->expectException(InvalidArgumentException::class);
        $rho->expectationValue($oper);
    }

    public function testStatevectorHelpers(): void
    {
        $state = $this->bellStatevector();

        $this->assertEqualsWithDelta(1.0, $state->densityMatrix()->purity(), 1.0e-12);
        $this->assertComplexMatrix([
            [0.5, 0.0],
            [0.0, 0.5],
        ], $state->reducedDensityMatrix([0])->data(), 1.0e-12);
        $this->assertComplexMatrix(
            $state->reducedDensityMatrix([0])->data(),
            $state->reduced_density_matrix([0])->data(),
            1.0e-12
        );
    }

    public function testStatevectorReducedDensityMatrixOrdering(): void
    {
        // 位相の異なる状態で、BLAS 経由 (Statevector) と
        // 部分トレース経由 (DensityMatrix) の縮約が一致するか確認する。
        $state = new Statevector(
            (new QuantumCircuit(2))->h(0)->cx(0, 1)->t(0)->s(1)
        );
        $rho = new DensityMatrix($state);

        foreach ([[0], [1], [0, 1], [1, 0]] as $keep) {
            $this->assertComplexMatrix(
                $rho->reducedDensityMatrix($keep)->data(),
                $state->reducedDensityMatrix($keep)->data(),
                1.0e-12
            );
        }
    }

    public function testThreeQubitReducedDensityMatrixOrdering(): void
    {
        // 3 qubit でも pure-state 経路 (transpose + gemm) と
        // 部分トレース経路 (transpose + gemv) が一致するか確認する。
        $state = new Statevector(
            (new QuantumCircuit(3))->h(0)->cx(0, 1)->t(1)->cx(1, 2)->s(0)->h(2)
        );
        $rho = new DensityMatrix($state);

        $orders = [
            [0], [1], [2],
            [0, 1], [1, 0], [0, 2], [2, 0], [1, 2], [2, 1],
            [0, 1, 2], [2, 1, 0], [1, 2, 0],
        ];
        foreach ($orders as $keep) {
            $this->assertComplexMatrix(
                $rho->reducedDensityMatrix($keep)->data(),
                $state->reducedDensityMatrix($keep)->data(),
                1.0e-12,
                'keep=' . implode(',', $keep)
            );
        }
        // partialTrace は残る qubit を元の昇順で並べる。
        $traceOrders = [[0], [1], [2], [0, 1], [0, 2], [1, 2], [0, 1, 2]];
        foreach ($traceOrders as $trace) {
            $keep = array_values(array_diff(range(0, 2), $trace));
            $this->assertComplexMatrix(
                $rho->partialTrace($trace)->data(),
                $state->reducedDensityMatrix($keep)->data(),
                1.0e-12,
                'trace=' . implode(',', $trace)
            );
        }
    }

    public function testDefaultDtypeFollowsPrecision(): void
    {
        $la = LinearAlgebraProvider::get();
        $rho = new DensityMatrix($la->array([
            C(1.0, 0.0),
            C(0.0, 0.0),
        ], dtype: NDArray::complex128));
        $this->assertSame(NDArray::complex128, $rho->ndarray()->dtype());

        try {
            Precision::setFloat32();
            $rho32 = new DensityMatrix([C(1.0, 0.0), C(0.0, 0.0)]);
            $this->assertSame(NDArray::complex64, $rho32->ndarray()->dtype());
            $this->assertSame(NDArray::complex64, DensityMatrix::fromInt(0, 1)->ndarray()->dtype());

            $reduced = $rho32->reducedDensityMatrix([0]);
            $this->assertSame(NDArray::complex64, $reduced->ndarray()->dtype());
        } finally {
            Precision::reset();
        }
    }

    public function testInvalidShapes(): void
    {
        $la = LinearAlgebraProvider::get();

        $this->expectException(InvalidArgumentException::class);
        new DensityMatrix($la->zeros([3, 3], NDArray::complex128));
    }

    public function testNonSquareMatrix(): void
    {
        $la = LinearAlgebraProvider::get();

        $this->expectException(InvalidArgumentException::class);
        new DensityMatrix($la->zeros([2, 3], NDArray::complex128));
    }

    public function testNonPowerOfTwoStatevector(): void
    {
        $la = LinearAlgebraProvider::get();

        $this->expectException(InvalidArgumentException::class);
        DensityMatrix::fromStatevector($la->zeros([3], NDArray::complex128));
    }

    public function testDuplicateQargs(): void
    {
        $rho = DensityMatrix::fromInt(0, 2);

        $this->expectException(InvalidArgumentException::class);
        $rho->partialTrace([0, 0]);
    }

    public function testOutOfRangeQargs(): void
    {
        $rho = DensityMatrix::fromInt(0, 2);

        $this->expectException(InvalidArgumentException::class);
        $rho->reducedDensityMatrix([2]);
    }

    public function testEmptyEnsemble(): void
    {
        $this->expectException(InvalidArgumentException::class);
        DensityMatrix::fromEnsemble([]);
    }

    /**
     * @param array<int,array<int,mixed>> $expected
     * @param array<int,array<int,mixed>> $actual
     */
    private function assertComplexMatrix(array $expected, array $actual, float $delta = 1.0e-12, string $message = ''): void
    {
        $prefix = $message === '' ? '' : $message . ': ';
        $this->assertCount(count($expected), $actual);
        foreach ($expected as $i => $row) {
            $this->assertCount(count($row), $actual[$i]);
            foreach ($row as $j => $value) {
                [$real, $imag] = $this->complexParts($value);
                $this->assertEqualsWithDelta(
                    $real,
                    (float)$actual[$i][$j]->real,
                    $delta,
                    "{$prefix}Real mismatch at ({$i},{$j})"
                );
                $this->assertEqualsWithDelta(
                    $imag,
                    (float)$actual[$i][$j]->imag,
                    $delta,
                    "{$prefix}Imag mismatch at ({$i},{$j})"
                );
            }
        }
    }

    /** @return array{float,float} */
    private function complexParts(mixed $value): array
    {
        if (is_int($value) || is_float($value)) {
            return [(float)$value, 0.0];
        }
        if (is_array($value)) {
            return [(float)($value[0] ?? 0.0), (float)($value[1] ?? 0.0)];
        }
        return [(float)($value->real ?? 0.0), (float)($value->imag ?? 0.0)];
    }

    /** @param array<int,float> $expected @param array<int,float> $actual */
    private function assertFloatVectorEquals(array $expected, array $actual, float $delta = 1.0e-6): void
    {
        $this->assertCount(count($expected), $actual);
        foreach ($expected as $i => $value) {
            $this->assertEqualsWithDelta($value, $actual[$i], $delta, "Mismatch at index {$i}.");
        }
    }
}
