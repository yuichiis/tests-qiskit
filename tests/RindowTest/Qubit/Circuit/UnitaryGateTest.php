<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Circuit;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Circuit\UnitaryGate;
use Rindow\Qubit\Math\LinearAlgebraProvider;
use Rindow\Math\Matrix\Complex;
use RindowTest\Qubit\TestUtils;

class UnitaryGateTest extends TestCase
{
    use TestUtils;

    public function testFromNumericArray(): void
    {
        $gate = new UnitaryGate([[0, 1], [1, 0]]);

        $this->assertSame('unitary', $gate->name());
        $this->assertSame(2, $gate->dimension());
        $this->assertGateMatrix($gate, [[0, 1], [1, 0]]);
    }

    public function testFromComplexValues(): void
    {
        $gate = new UnitaryGate([
            [new Complex(0.0, 0.0), new Complex(0.0, -1.0)],
            [['real' => 0.0, 'imag' => 1.0], (object)['real' => 1.0, 'imag' => 0.0]],
        ]);

        $this->assertSame(2, $gate->dimension());
        $this->assertGateMatrix($gate, [
            [0, ['real' => 0.0, 'imag' => -1.0]],
            [['real' => 0.0, 'imag' => 1.0], 1],
        ]);
    }

    #[DataProvider('ndarrayDtypeProvider')]
    public function testFromNDArray(int $dtype): void
    {
        $la = LinearAlgebraProvider::get();
        $matrix = $la->array([[1, 0], [0, -1]], dtype: $dtype);
        $gate = new UnitaryGate($matrix);

        $this->assertSame(2, $gate->dimension());
        $this->assertGateMatrix($gate, [[1, 0], [0, -1]]);
    }

    /** @return array<string,array{int}> */
    public static function ndarrayDtypeProvider(): array
    {
        return [
            'int32' => [NDArray::int32],
            'float32' => [NDArray::float32],
            'float64' => [NDArray::float64],
            'complex128' => [NDArray::complex128],
        ];
    }

    public function testFourDimensionalMatrix(): void
    {
        $gate = new UnitaryGate($this->cnotMatrix());

        $this->assertSame(4, $gate->dimension());
        $this->assertGateMatrix($gate, $this->cnotMatrix());
    }

    public function testControlProducesControlledGate(): void
    {
        $controlled = (new UnitaryGate([[0, 1], [1, 0]]))->control();

        $this->assertInstanceOf(UnitaryGate::class, $controlled);
        $this->assertSame(4, $controlled->dimension());
        $this->assertGateMatrix($controlled, $this->cnotMatrix());
    }

    public function testControlRequiresPositiveCount(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new UnitaryGate([[0, 1], [1, 0]]))->control(0);
    }

    public function testFredkinGateMatrix(): void
    {
        $gate = UnitaryGate::fredkin();
        $expected = [];
        for ($i = 0; $i < 8; $i++) {
            $expected[$i] = array_fill(0, 8, 0);
        }
        $map = [0 => 0, 1 => 1, 2 => 2, 3 => 5, 4 => 4, 5 => 3, 6 => 6, 7 => 7];
        foreach ($map as $input => $output) {
            $expected[$output][$input] = 1;
        }

        $this->assertSame(8, $gate->dimension());
        $this->assertGateMatrix($gate, $expected);
    }

    public function testDraw(): void
    {
        $this->assertSame(
            ['cells' => [0 => 'unitary'], 'span' => [0]],
            (new UnitaryGate([[0, 1], [1, 0]]))->draw([0])
        );
        $this->assertSame(
            ['cells' => [0 => 'U', 1 => 'U'], 'span' => [0, 1]],
            (new UnitaryGate($this->cnotMatrix()))->draw([0, 1])
        );
        $this->assertEquals(
            ['cells' => [0 => 'U', 2 => 'U', 1 => '│'], 'span' => [0, 1, 2]],
            (new UnitaryGate($this->cnotMatrix()))->draw([0, 2])
        );
    }

    #[DataProvider('invalidMatrixProvider')]
    public function testInvalidMatrix(array|NDArray $matrix): void
    {
        $this->expectException(InvalidArgumentException::class);
        new UnitaryGate($matrix);
    }

    /** @return array<string,array{array|NDArray}> */
    public static function invalidMatrixProvider(): array
    {
        $la = LinearAlgebraProvider::get();
        return [
            'non square array' => [[[1, 0, 1], [0, 1, 0]]],
            'one by one' => [[[1]]],
            'three by three' => [[
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
            ]],
            'non square NDArray' => [$la->array([[1, 0, 0], [0, 1, 0]], dtype: NDArray::complex64)],
            'plain re/im pair' => [[[[0, 1], [1, 0]], [[1, 0], [0, 1]]]],
            'unsupported object' => [[[new \stdClass(), new \stdClass()], [new \stdClass(), new \stdClass()]]],
        ];
    }

    /** @return array<array<int>> */
    private function cnotMatrix(): array
    {
        return [
            [1, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0],
            [0, 1, 0, 0],
        ];
    }
}
