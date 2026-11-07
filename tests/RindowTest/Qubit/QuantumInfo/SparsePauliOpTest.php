<?php
declare(strict_types=1);

namespace RindowTest\Qubit\QuantumInfo;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Math\LinearAlgebraProvider;
use Rindow\Qubit\QuantumInfo\Pauli;
use Rindow\Qubit\QuantumInfo\SparsePauliOp;
use Rindow\Qubit\QuantumInfo\Statevector;
use Rindow\Qubit\QuantumCircuit;
use function Rindow\Math\Matrix\C;

class SparsePauliOpTest extends TestCase
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

    public function testPauliLabelAndQubitOrder(): void
    {
        $pauli = new Pauli('IZ');

        $this->assertSame('IZ', $pauli->label());
        $this->assertSame(2, $pauli->numQubits());
        $this->assertSame('I', $pauli->pauliAt(1));
        $this->assertSame('Z', $pauli->pauliAt(0));
    }

    public function testInvalidPauliLabel(): void
    {
        $this->expectException(InvalidArgumentException::class);
        new Pauli('IA');
    }

    public function testSingleQubitPauliMatrices(): void
    {
        $i = [[1, 0], [0, 1]];
        $x = [[0, 1], [1, 0]];
        $y = [[0, [0.0, -1.0]], [[0.0, 1.0], 0]];
        $z = [[1, 0], [0, -1]];

        $this->assertComplexMatrix($i, (new Pauli('I'))->toMatrix()->toArray());
        $this->assertComplexMatrix($x, (new Pauli('X'))->toMatrix()->toArray());
        $this->assertComplexMatrix($y, (new Pauli('Y'))->toMatrix()->toArray());
        $this->assertComplexMatrix($z, (new Pauli('Z'))->toMatrix()->toArray());
    }

    public function testTwoQubitPauliMatrices(): void
    {
        // IZ = Z on qubit 0 (LSB), ZI = Z on qubit 1
        $this->assertComplexMatrix(
            [[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]],
            (new Pauli('IZ'))->toMatrix()->toArray()
        );
        $this->assertComplexMatrix(
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]],
            (new Pauli('ZI'))->toMatrix()->toArray()
        );
        $this->assertComplexMatrix(
            [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]],
            (new Pauli('XX'))->toMatrix()->toArray()
        );
    }

    public function testFromListAndToList(): void
    {
        $op = SparsePauliOp::fromList([['II', 1.5], ['ZZ', -0.5]]);

        $this->assertSame(2, $op->numQubits());
        $this->assertSame(2, $op->num_qubits());
        $this->assertSame(2, $op->numTerms());
        $this->assertSame([['II', 1.5], ['ZZ', -0.5]], $op->toList());
        $this->assertSame($op->toList(), $op->to_list());
        $this->assertSame($op->terms(), $op->terms());
    }

    public function testLabelLengthMismatch(): void
    {
        $this->expectException(InvalidArgumentException::class);
        SparsePauliOp::fromList([['II', 1.0], ['ZZZ', 1.0]]);
    }

    public function testH2MatrixMatchesQiskit(): void
    {
        $matrix64 = $this->h2()->toMatrix();
        $matrix = $matrix64->toArray();

        $this->assertSame(\Interop\Polite\Math\Matrix\NDArray::complex128, $matrix64->dtype());
        $this->assertSame([4, 4], $matrix64->shape());
        $this->assertComplexMatrix([
            [-1.0636533500290943, 0.0, 0.0, 0.18093119978423156],
            [0.0, -1.8369679912029844, 0.18093119978423156, 0.0],
            [0.0, 0.18093119978423156, -0.24521829183026272, 0.0],
            [0.18093119978423156, 0.0, 0.0, -1.0636533500290943],
        ], $matrix);
    }

    public function testToMatrixAddsAndSubtractsNegativeCoefficients(): void
    {
        // 正の係数は加算、負の係数は減算として axpy で積み上げられる。
        $op = SparsePauliOp::fromList([
            ['IZ', 1.5],
            ['ZI', -0.25],
            ['XX', 0.5],
        ]);
        $expected = [
            [1.25, 0.0, 0.0, 0.5],
            [0.0, -1.75, 0.5, 0.0],
            [0.0, 0.5, 1.75, 0.0],
            [0.5, 0.0, 0.0, -1.25],
        ];

        $this->assertComplexMatrix($expected, $op->toMatrix()->toArray(), 1.0e-6);
        $this->assertComplexMatrix(
            $expected,
            $op->toMatrix(\Interop\Polite\Math\Matrix\NDArray::complex128)->toArray()
        );
    }

    public function testPackedExpectationValueMatchesDenseMatrix(): void
    {
        // Y を含む項は複素位相を持つ。上三角パックド行列 + hpmv が
        // 密行列 + gemv と同じ期待値を返すことを確認する。
        $op = SparsePauliOp::fromList([
            ['II', -0.5],
            ['XI', 0.35],
            ['YI', -0.2],
            ['IY', 0.15],
            ['XX', 0.18],
            ['YY', 0.11],
            ['ZZ', -0.07],
        ]);
        $la = LinearAlgebraProvider::get();
        $amplitudes = [C(0.5, 0.1), C(-0.2, 0.3), C(0.4, -0.25), C(0.0, 0.6)];

        foreach ([NDArray::complex64, NDArray::complex128] as $dtype) {
            $matrix = $op->toMatrix($dtype);
            $state = $la->array($amplitudes, dtype: $dtype);
            $expected = $la->scalar($la->dot($state, $la->gemv($matrix, $state)));
            $actual = $op->expectationValueNDArray($state);

            $delta = $dtype === NDArray::complex64 ? 1.0e-5 : 1.0e-12;
            $this->assertEqualsWithDelta((float)$expected->real, $actual, $delta);
            // エルミート演算子なので虚部は数値誤差だけ。
            $this->assertEqualsWithDelta(0.0, (float)$expected->imag, $delta);
        }
    }

    public function testPackedExpectationValueForThreeQubits(): void
    {
        // 4x4 を超えるサイズでもパックド位置が正しいことを確認する。
        $op = SparsePauliOp::fromList([
            ['IIX', 0.3],
            ['IYI', -0.4],
            ['ZIZ', 0.2],
            ['XXX', 0.1],
            ['YIY', -0.05],
        ]);
        $la = LinearAlgebraProvider::get();
        $amplitudes = [];
        for ($i = 0; $i < 8; $i++) {
            $amplitudes[] = C(0.1 * $i - 0.4, 0.05 * (($i * 7) % 8) - 0.2);
        }
        $state = $la->array($amplitudes, dtype: NDArray::complex128);
        $matrix = $op->toMatrix();
        $expected = $la->scalar($la->dot($state, $la->gemv($matrix, $state)));

        $this->assertEqualsWithDelta((float)$expected->real, $op->expectationValueNDArray($state), 1.0e-12);
        $this->assertEqualsWithDelta(0.0, (float)$expected->imag, 1.0e-12);
    }

    public function testExpectationValueOnComputationalBasis(): void
    {
        $h = $this->h2();

        $this->assertEqualsWithDelta(
            -1.0636533500290943,
            $h->expectationValue(new Statevector($this->basisState(2, 0))),
            1.0e-6
        );
        $this->assertEqualsWithDelta(
            -1.8369679912029844,
            $h->expectationValue(new Statevector($this->basisState(2, 1))),
            1.0e-6
        );
    }

    public function testExpectationValueMatchesMatrix(): void
    {
        $h = $this->h2();
        $circuit = (new QuantumCircuit(2))->h(0)->cx(0, 1);
        $state = new Statevector($circuit);

        $matrix = $h->toMatrix(\Interop\Polite\Math\Matrix\NDArray::complex128)->toArray();
        $amplitudes = $state->data();
        $expectedRe = 0.0;
        $expectedIm = 0.0;
        for ($row = 0; $row < 4; $row++) {
            for ($col = 0; $col < 4; $col++) {
                $a = $amplitudes[$row];
                $b = $amplitudes[$col];
                $m = $matrix[$row][$col];
                $expectedRe += $a->real * $m->real * $b->real
                    - $a->real * $m->imag * $b->imag
                    + $a->imag * $m->real * $b->imag
                    + $a->imag * $m->imag * $b->real;
                $expectedIm += $a->real * $m->real * $b->imag
                    + $a->real * $m->imag * $b->real
                    - $a->imag * $m->real * $b->real
                    + $a->imag * $m->imag * $b->imag;
            }
        }

        // 状態ベクトルと演算子行列は同じ既定精度 (complex128) を使う。
        $this->assertEqualsWithDelta($expectedRe, $h->expectationValue($state), 1.0e-12);
        $this->assertEqualsWithDelta(0.0, $expectedIm, 1.0e-12);
        $this->assertEqualsWithDelta($expectedRe, $h->expectation_value($circuit), 1.0e-12);
    }

    public function testStatevectorQubitCountMismatch(): void
    {
        $this->expectException(InvalidArgumentException::class);
        $this->h2()->expectationValue(new Statevector($this->basisState(1, 0)));
    }

    public function testExpectationValueAcceptsNDArrayDirectly(): void
    {
        $h = $this->h2();
        $sim = new \Rindow\Qubit\Quantum\StatevectorSimulator(2);
        $sim->applyRy(0, 0.3);
        $sim->applyRy(1, -0.4);
        $state = $sim->stateNDArray();

        $this->assertSame(
            $state->dtype(),
            $h->toMatrix()->dtype(),
            'The BLAS path should keep the state and matrix dtypes aligned.'
        );
        $this->assertEqualsWithDelta(
            $h->expectationValue(new Statevector($state)),
            $h->expectationValueNDArray($state),
            1.0e-12
        );
    }

    public function testToMatrixReturnsIndependentCopies(): void
    {
        $h = $this->h2();
        $matrix = $h->toMatrix();
        $dtype = $matrix->dtype();

        // 返り値を書き換えても内部キャッシュ (期待値計算) に影響しないこと。
        $matrix[0][0] = \Rindow\Math\Matrix\C(999.0, 0.0);

        $again = $h->toMatrix();
        $this->assertSame($dtype, $again->dtype());
        $this->assertEqualsWithDelta(-1.0636533500290943, $again[0][0]->real, 1.0e-6);

        $energy = $h->expectationValue(new Statevector($this->basisState(2, 0)));
        $this->assertEqualsWithDelta(-1.0636533500290943, $energy, 1.0e-6);
    }

    public function testMatrixCacheIsReusedAcrossCalls(): void
    {
        $h = $this->h2();

        $matrix64 = $h->toMatrix(\Interop\Polite\Math\Matrix\NDArray::complex64);
        $this->assertSame(
            \Interop\Polite\Math\Matrix\NDArray::complex64,
            $matrix64->dtype()
        );
        $this->assertSame(
            \Interop\Polite\Math\Matrix\NDArray::complex128,
            $h->toMatrix()->dtype()
        );
        // dtype ごとにキャッシュが分かれていること。
        $this->assertSame(
            \Interop\Polite\Math\Matrix\NDArray::complex64,
            $h->toMatrix(\Interop\Polite\Math\Matrix\NDArray::complex64)->dtype()
        );
    }

    /**
     * @return array<int,object{real:float,imag:float}>
     */
    private function basisState(int $numQubits, int $index): array
    {
        $state = [];
        for ($i = 0; $i < (1 << $numQubits); $i++) {
            $state[] = (object)[
                'real' => $i === $index ? 1.0 : 0.0,
                'imag' => 0.0,
            ];
        }
        return $state;
    }

    /**
     * @param array<int,array<int,mixed>> $expected
     * @param array<int,array<int,object{real:float,imag:float}>> $actual
     */
    private function assertComplexMatrix(array $expected, array $actual, float $delta = 1.0e-12): void
    {
        $this->assertCount(count($expected), $actual);
        foreach ($expected as $i => $row) {
            $this->assertCount(count($row), $actual[$i]);
            foreach ($row as $j => $value) {
                [$real, $imag] = $this->complexParts($value);
                $this->assertEqualsWithDelta($real, $actual[$i][$j]->real, $delta, "Real mismatch at ({$i},{$j})");
                $this->assertEqualsWithDelta($imag, $actual[$i][$j]->imag, $delta, "Imag mismatch at ({$i},{$j})");
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
}
