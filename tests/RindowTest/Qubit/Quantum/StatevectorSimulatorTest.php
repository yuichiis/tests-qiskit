<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Quantum;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use PHPUnit\Framework\TestCase;
use Rindow\Math\Matrix\Complex;
use Rindow\Qubit\Circuit\QFTGate;
use Rindow\Qubit\Circuit\UnitaryGate;
use Rindow\Qubit\Quantum\StatevectorSimulator;
use Rindow\Qubit\QuantumCircuit;
use RindowTest\Qubit\TestUtils;

class StatevectorSimulatorTest extends TestCase
{
    use TestUtils;

    public function testConstructorAndInitialState(): void
    {
        $sim = new StatevectorSimulator(2);

        $this->assertSame(2, $sim->numQubits());
        $this->assertSame(4, $sim->dimension());
        $this->assertAmplitudesEqual([[1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]], $sim->getStatevector());
        $this->assertInstanceOf(NDArray::class, $sim->stateNDArray());
        $this->assertSame([4], $sim->stateNDArray()->shape());
        $this->assertNotNull($sim->la());
    }

    public function testInvalidQubitCount(): void
    {
        $this->expectException(InvalidArgumentException::class);
        new StatevectorSimulator(0);
    }

    public function testTooManyQubits(): void
    {
        $this->expectException(InvalidArgumentException::class);
        new StatevectorSimulator(31);
    }

    public function testApplyCircuitQubitMismatch(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new StatevectorSimulator(1))->applyCircuit(new QuantumCircuit(2));
    }

    public function testHardamardAndPauliGates(): void
    {
        $sim = new StatevectorSimulator(1);
        $sim->applyOperation('h', [0]);
        $this->assertAmplitudesEqual(
            [[1.0 / sqrt(2.0), 0.0], [1.0 / sqrt(2.0), 0.0]],
            $sim->getStatevector()
        );

        $sim = new StatevectorSimulator(1);
        $sim->applyOperation('x', [0]);
        $this->assertAmplitudesEqual([[0.0, 0.0], [1.0, 0.0]], $sim->getStatevector());

        $sim = new StatevectorSimulator(1);
        $sim->applyOperation('y', [0]);
        $this->assertAmplitudesEqual([[0.0, 0.0], [0.0, 1.0]], $sim->getStatevector());

        $sim = new StatevectorSimulator(1);
        $sim->applyOperation('z', [0]);
        $this->assertAmplitudesEqual([[1.0, 0.0], [0.0, 0.0]], $sim->getStatevector());
    }

    public function testBellState(): void
    {
        $qc = (new QuantumCircuit(2))->h(0)->cx(0, 1);
        $sim = new StatevectorSimulator(2);
        $sim->applyCircuit($qc);

        $this->assertAmplitudesEqual(
            [[1.0 / sqrt(2.0), 0.0], [0.0, 0.0], [0.0, 0.0], [1.0 / sqrt(2.0), 0.0]],
            $sim->getStatevector()
        );
        $this->assertFloatVectorEquals(
            [0.5, 0.0, 0.0, 0.5],
            $sim->getProbabilities()->toArray()
        );
    }

    public function testEndiannessOfQubits(): void
    {
        $sim = new StatevectorSimulator(2);
        $sim->applyOperation('x', [1]);
        $this->assertAmplitudesEqual(
            [[0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 0.0]],
            $sim->getStatevector()
        );
    }

    public function testCxControlOrdering(): void
    {
        $sim = new StatevectorSimulator(2);
        $sim->applyOperation('x', [1]);
        $sim->applyOperation('cx', [0, 1]);
        $this->assertAmplitudesEqual(
            [[0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 0.0]],
            $sim->getStatevector()
        );

        $sim = new StatevectorSimulator(2);
        $sim->applyOperation('x', [1]);
        $sim->applyOperation('cx', [1, 0]);
        $this->assertAmplitudesEqual(
            [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [1.0, 0.0]],
            $sim->getStatevector()
        );

        $sim = new StatevectorSimulator(2);
        $sim->applyOperation('x', [0]);
        $sim->applyOperation('cx', [0, 1]);
        $this->assertAmplitudesEqual(
            [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [1.0, 0.0]],
            $sim->getStatevector()
        );
    }

    public function testCzAndCp(): void
    {
        $sim = new StatevectorSimulator(2);
        $sim->applyOperation('x', [0]);
        $sim->applyOperation('x', [1]);
        $sim->applyOperation('cz', [0, 1]);
        $this->assertAmplitudesEqual(
            [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [-1.0, 0.0]],
            $sim->getStatevector()
        );

        $sim = new StatevectorSimulator(2);
        $sim->applyOperation('x', [0]);
        $sim->applyOperation('x', [1]);
        $sim->applyOperation('cp', [0, 1], [M_PI / 2.0]);
        $this->assertAmplitudesEqual(
            [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 1.0]],
            $sim->getStatevector()
        );
    }

    public function testSwap(): void
    {
        $sim = new StatevectorSimulator(2);
        $sim->applyOperation('x', [0]);
        $sim->applyOperation('swap', [0, 1]);
        $this->assertAmplitudesEqual(
            [[0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 0.0]],
            $sim->getStatevector()
        );
    }

    public function testRotationAndPhaseHelpers(): void
    {
        $sim = new StatevectorSimulator(1);
        $sim->applyRx(0, M_PI / 2.0);
        $this->assertAmplitudesEqual(
            [[1.0 / sqrt(2.0), 0.0], [0.0, -1.0 / sqrt(2.0)]],
            $sim->getStatevector()
        );

        $sim = new StatevectorSimulator(1);
        $sim->applyRy(0, M_PI / 2.0);
        $this->assertAmplitudesEqual(
            [[1.0 / sqrt(2.0), 0.0], [1.0 / sqrt(2.0), 0.0]],
            $sim->getStatevector()
        );

        $sim = new StatevectorSimulator(1);
        $sim->applyRz(0, M_PI / 2.0);
        $this->assertAmplitudesEqual(
            [[1.0 / sqrt(2.0), -1.0 / sqrt(2.0)], [0.0, 0.0]],
            $sim->getStatevector()
        );

        $sim = new StatevectorSimulator(1);
        $sim->applyOperation('x', [0]);
        $sim->applyPhase(M_PI / 2.0, 0);
        $this->assertAmplitudesEqual([[0.0, 0.0], [0.0, 1.0]], $sim->getStatevector());

        $sim = new StatevectorSimulator(1);
        $sim->applyU(M_PI / 2.0, M_PI / 3.0, M_PI / 4.0, 0);
        $this->assertAmplitudesEqual(
            [
                [1.0 / sqrt(2.0), 0.0],
                [0.5 * sqrt(0.5), 0.5 * sqrt(1.5)],
            ],
            $sim->getStatevector()
        );
    }

    public function testApplyNQubitOnArbitraryQubits(): void
    {
        $swap = [
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 1, 0, 0],
            [0, 0, 0, 1],
        ];

        $sim = new StatevectorSimulator(3);
        $sim->applyOperation('x', [0]);
        // First argument is MSB: swap q2 and q0.
        $sim->applyNQubit([2, 0], $swap);

        $this->assertAmplitudesEqual(
            [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
            $sim->getStatevector()
        );
    }

    public function testFredkinGateThroughCircuit(): void
    {
        $qc = new QuantumCircuit(3);
        $qc->x(0)->x(1);
        $qc->append(UnitaryGate::fredkin(), [0, 1, 2]);

        $sim = new StatevectorSimulator(3);
        $sim->applyCircuit($qc);

        $probabilities = $sim->getProbabilities()->toArray();
        $this->assertEqualsWithDelta(1.0, $probabilities[5], 1.0e-6);
        for ($i = 0; $i < 8; $i++) {
            if ($i !== 5) {
                $this->assertEqualsWithDelta(0.0, $probabilities[$i], 1.0e-6);
            }
        }
    }

    public function testQftAndInverseRestoreState(): void
    {
        $qft = new QFTGate(2);
        $qc = new QuantumCircuit(2);
        $qc->x(0);
        $qc->append($qft, [0, 1]);
        $qc->append($qft->inverse(), [0, 1]);

        $sim = new StatevectorSimulator(2);
        $sim->applyCircuit($qc);

        $this->assertAmplitudesEqual(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
            $sim->getStatevector(),
            1.0e-5
        );
    }

    public function testApplyOperationUnitaryLegacyFormat(): void
    {
        $sim = new StatevectorSimulator(1);
        $sim->applyOperation('unitary', [0], [
            new Complex(0.0, 0.0),
            new Complex(1.0, 0.0),
            new Complex(1.0, 0.0),
            new Complex(0.0, 0.0),
        ]);

        $this->assertAmplitudesEqual([[0.0, 0.0], [1.0, 0.0]], $sim->getStatevector());
    }

    public function testApplyOperationUnitaryNDArray(): void
    {
        $la = \Rindow\Qubit\Math\LinearAlgebraProvider::get();
        $matrix = $la->array([[0, 1], [1, 0]], dtype: NDArray::complex64);

        $sim = new StatevectorSimulator(1);
        $sim->applyOperation('unitary', [0], $matrix);

        $this->assertAmplitudesEqual([[0.0, 0.0], [1.0, 0.0]], $sim->getStatevector());
    }

    public function testApplyOperationUnitaryNonComplexNDArray(): void
    {
        $la = \Rindow\Qubit\Math\LinearAlgebraProvider::get();
        $matrix = $la->array([[0, 1], [1, 0]], dtype: NDArray::float64);

        $sim = new StatevectorSimulator(1);
        $sim->applyOperation('unitary', [0], $matrix);

        $this->assertAmplitudesEqual([[0.0, 0.0], [1.0, 0.0]], $sim->getStatevector());
    }

    public function testMeasureBarrierAndResetAreNoOps(): void
    {
        $sim = new StatevectorSimulator(1);
        $sim->applyOperation('measure', [0]);
        $sim->applyOperation('barrier', [0]);
        $sim->applyOperation('reset', [0]);

        $this->assertAmplitudesEqual([[1.0, 0.0], [0.0, 0.0]], $sim->getStatevector());
    }

    public function testUnknownGate(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new StatevectorSimulator(1))->applyOperation('unknown', [0]);
    }

    public function testApplySingleValidation(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new StatevectorSimulator(1))->applySingle(1, [[1, 0], [0, 1]]);
    }

    public function testApplySingleMatrixSizeValidation(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new StatevectorSimulator(2))->applySingle(0, [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]);
    }

    public function testApplyNQubitValidation(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new StatevectorSimulator(1))->applyNQubit([], [[1]]);
    }

    public function testApplyNQubitMatrixSizeValidation(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new StatevectorSimulator(2))->applyNQubit([1, 0], [[1, 0], [0, 1]]);
    }

    public function testApplyNQubitDuplicateQubitValidation(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new StatevectorSimulator(2))->applyNQubit([0, 0], [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]);
    }

    public function testProbabilitiesWithMeasurementMap(): void
    {
        $sim = new StatevectorSimulator(2);
        $sim->applyOperation('x', [0]);

        $this->assertSame([1 => 1.0], $sim->probabilities([0 => 0], 1));
        $this->assertSame([10 => 1.0], $sim->probabilities([0 => 1], 2));
        $this->assertSame(['01' => 1.0], $sim->probabilities());
    }

    public function testSampleCountsExactState(): void
    {
        $sim = new StatevectorSimulator(1);
        $sim->applyOperation('x', [0]);

        $this->assertSame([1 => 100], $sim->sampleCounts([0 => 0], 1, 100, 42));
    }

    public function testSampleCountsIsSeededAndSummed(): void
    {
        $sim = new StatevectorSimulator(2);
        $sim->applyOperation('h', [0]);
        $sim->applyOperation('cx', [0, 1]);

        $first = $sim->sampleCounts([0 => 0, 1 => 1], 2, 200, 1234);
        $second = $sim->sampleCounts([0 => 0, 1 => 1], 2, 200, 1234);

        $this->assertSame($first, $second);
        $this->assertSame(200, array_sum($first));
        foreach (array_keys($first) as $bits) {
            $this->assertContains((string)$bits, ['00', '11']);
        }
    }

    public function testSampleCountsWithMeasurementMap(): void
    {
        $sim = new StatevectorSimulator(1);
        $sim->applyOperation('x', [0]);

        $this->assertSame([10 => 25], $sim->sampleCounts([0 => 1], 2, 25, 7));
    }

    public function testSampleCountsRequiresShots(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new StatevectorSimulator(1))->sampleCounts([0 => 0], 1, 0);
    }
}
