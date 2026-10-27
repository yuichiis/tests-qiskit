<?php
declare(strict_types=1);

namespace RindowTest\Qubit\QuantumInfo;

use InvalidArgumentException;
use PHPUnit\Framework\TestCase;
use Rindow\Math\Matrix\Complex;
use Rindow\Qubit\Quantum\StatevectorSimulator;
use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\QuantumInfo\Statevector;
use RindowTest\Qubit\TestUtils;

class StatevectorTest extends TestCase
{
    use TestUtils;

    public function testConstructFromArray(): void
    {
        $statevector = new Statevector([
            new Complex(1.0, 0.0),
            new Complex(0.0, 0.0),
            new Complex(0.0, 0.0),
            new Complex(0.0, 0.0),
        ]);

        $this->assertSame(2, $statevector->numQubits());
        $this->assertSame(2, $statevector->num_qubits());
        $this->assertSame(4, $statevector->dim());
        $this->assertSame(4, $statevector->dimension());
        $this->assertCount(4, $statevector);
        $this->assertCount(4, $statevector->data());
        $this->assertSame($statevector->data(), $statevector->getStatevector());
        $this->assertAmplitudesEqual(
            [[1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
            $statevector->data()
        );
    }

    public function testConstructFromPlainObjects(): void
    {
        $statevector = new Statevector([
            (object)['real' => 0.0, 'imag' => 0.0],
            (object)['real' => 1.0, 'imag' => 0.0],
        ]);

        $this->assertSame(1, $statevector->numQubits());
        $this->assertSame(2, $statevector->dim());
        $this->assertAmplitudesEqual([[0.0, 0.0], [1.0, 0.0]], $statevector->data());
    }

    public function testConstructFromArrayNormalizesKeys(): void
    {
        $statevector = new Statevector([
            2 => (object)['real' => 1.0, 'imag' => 0.0],
            0 => (object)['real' => 0.0, 'imag' => 0.0],
            1 => (object)['real' => 0.0, 'imag' => 0.0],
            3 => (object)['real' => 0.0, 'imag' => 0.0],
        ]);

        $this->assertAmplitudesEqual(
            [[1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
            $statevector->data()
        );
    }

    public function testConstructFromCircuit(): void
    {
        $qc = (new QuantumCircuit(2))->h(0)->cx(0, 1);
        $statevector = new Statevector($qc);

        $this->assertSame(2, $statevector->numQubits());
        $this->assertAmplitudesEqual(
            [[1.0 / sqrt(2.0), 0.0], [0.0, 0.0], [0.0, 0.0], [1.0 / sqrt(2.0), 0.0]],
            $statevector->data()
        );
    }

    public function testConstructFromSimulator(): void
    {
        $sim = new StatevectorSimulator(1);
        $sim->applyOperation('x', [0]);
        $statevector = new Statevector($sim);

        $this->assertSame(1, $statevector->numQubits());
        $this->assertSame(2, $statevector->dim());
        $this->assertAmplitudesEqual([[0.0, 0.0], [1.0, 0.0]], $statevector->data());
    }

    public function testCopyConstructor(): void
    {
        $original = new Statevector([
            new Complex(1.0, 0.0),
            new Complex(0.0, 0.0),
        ]);
        $copy = new Statevector($original);

        $this->assertSame($original->data(), $copy->data());
        $this->assertSame($original->numQubits(), $copy->numQubits());
        $this->assertSame($original->dim(), $copy->dim());
    }

    public function testIteration(): void
    {
        $statevector = new Statevector([
            new Complex(1.0, 0.0),
            new Complex(0.0, 0.0),
        ]);

        $values = [];
        foreach ($statevector as $value) {
            $values[] = $value;
        }

        $this->assertSame($statevector->data(), $values);
    }

    public function testInvalidLength(): void
    {
        $this->expectException(InvalidArgumentException::class);
        new Statevector([
            new Complex(1.0, 0.0),
            new Complex(0.0, 0.0),
            new Complex(0.0, 0.0),
        ]);
    }

    public function testEmptyLength(): void
    {
        $this->expectException(InvalidArgumentException::class);
        new Statevector([]);
    }

    public function testProbabilities(): void
    {
        $statevector = new Statevector([
            new Complex(1.0 / sqrt(2.0), 0.0),
            new Complex(0.0, 0.0),
            new Complex(0.0, 0.0),
            new Complex(0.0, 1.0 / sqrt(2.0)),
        ]);

        $this->assertFloatVectorEquals([0.5, 0.0, 0.0, 0.5], $statevector->probabilities());
        $this->assertFloatVectorEquals([0.5, 0.0, 0.0, 0.5], $statevector->probabilities(null, 3));
    }

    public function testPartialProbabilities(): void
    {
        $statevector = new Statevector([
            new Complex(0.0, 0.0),
            new Complex(1.0, 0.0),
            new Complex(0.0, 0.0),
            new Complex(0.0, 0.0),
        ]);

        $this->assertFloatVectorEquals([0.0, 1.0], $statevector->probabilities([0]));
        $this->assertFloatVectorEquals([1.0, 0.0], $statevector->probabilities([1]));
        $this->assertFloatVectorEquals([0.0, 0.0, 1.0, 0.0], $statevector->probabilities([1, 0]));
        $this->assertFloatVectorEquals([0.0, 1.0, 0.0, 0.0], $statevector->probabilities([0, 1]));
    }

    public function testTooManyQargs(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new Statevector([new Complex(1.0, 0.0), new Complex(0.0, 0.0)]))->probabilities([0, 1, 0]);
    }

    public function testOutOfRangeQarg(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new Statevector([new Complex(1.0, 0.0), new Complex(0.0, 0.0)]))->probabilities([2]);
    }

    public function testDuplicateQargs(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new Statevector([new Complex(1.0, 0.0), new Complex(0.0, 0.0)]))->probabilities([0, 0]);
    }

    public function testProbabilitiesDict(): void
    {
        $statevector = new Statevector([
            new Complex(1.0 / sqrt(2.0), 0.0),
            new Complex(0.0, 0.0),
            new Complex(0.0, 0.0),
            new Complex(1.0 / sqrt(2.0), 0.0),
        ]);

        $dict = $statevector->probabilitiesDict();
        $this->assertCount(2, $dict);
        $this->assertArrayHasKey('00', $dict);
        $this->assertArrayHasKey(11, $dict);
        $this->assertEqualsWithDelta(0.5, $dict['00'], 1.0e-6);
        $this->assertEqualsWithDelta(0.5, $dict[11], 1.0e-6);
        $this->assertSame($dict, $statevector->probabilities_dict());
    }

    public function testProbabilitiesDictExcludesZeroProbability(): void
    {
        $statevector = new Statevector([
            new Complex(1.0, 0.0),
            new Complex(0.0, 0.0),
            new Complex(0.0, 0.0),
            new Complex(0.0, 0.0),
        ]);

        $this->assertSame(['00' => 1.0], $statevector->probabilitiesDict());
    }

    public function testProbabilitiesDictWithQargs(): void
    {
        $statevector = new Statevector([
            new Complex(0.0, 0.0),
            new Complex(1.0, 0.0),
            new Complex(0.0, 0.0),
            new Complex(0.0, 0.0),
        ]);

        $this->assertSame([1 => 1.0], $statevector->probabilitiesDict([0]));
        $this->assertSame([0 => 1.0], $statevector->probabilitiesDict([1]));
    }
}
