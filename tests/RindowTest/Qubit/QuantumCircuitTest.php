<?php
declare(strict_types=1);

namespace RindowTest\Qubit;

use ArrayIterator;
use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Circuit\HGate;
use Rindow\Qubit\QuantumCircuit;

class QuantumCircuitTest extends TestCase
{
    public function testConstructorAndAccessors(): void
    {
        $qc = new QuantumCircuit(3);
        $this->assertSame(3, $qc->numQubits());
        $this->assertSame(0, $qc->numClbits());
        $this->assertSame(3, $qc->num_qubits());
        $this->assertSame(0, $qc->num_clbits());

        $qc = new QuantumCircuit(2, 3);
        $this->assertSame(2, $qc->numQubits());
        $this->assertSame(3, $qc->numClbits());
    }

    public function testInvalidConstructorArguments(): void
    {
        $this->expectException(InvalidArgumentException::class);
        new QuantumCircuit(0);
    }

    public function testTooManyQubits(): void
    {
        $this->expectException(InvalidArgumentException::class);
        new QuantumCircuit(31);
    }

    public function testNegativeClbits(): void
    {
        $this->expectException(InvalidArgumentException::class);
        new QuantumCircuit(1, -1);
    }

    public function testSingleQubitGateHelpers(): void
    {
        $qc = new QuantumCircuit(3);
        $qc->h([0, 1])
            ->x(2)
            ->y(new ArrayIterator([2, 0]))
            ->z(1)
            ->s(0)
            ->sdg(0)
            ->t(0)
            ->tdg(0)
            ->id(0);

        $this->assertSame(
            ['h', 'h', 'x', 'y', 'y', 'z', 's', 'sdg', 't', 'tdg', 'id'],
            array_column($qc->operations(), 'name')
        );
        $this->assertSame(
            [[0], [1], [2], [2], [0], [1], [0], [0], [0], [0], [0]],
            array_column($qc->operations(), 'qubits')
        );
    }

    public function testParametrizedSingleQubitGateHelpers(): void
    {
        $qc = new QuantumCircuit(1);
        $qc->rx(0.1, 0)
            ->ry(0.2, 0)
            ->rz(0.3, 0)
            ->p(0.4, 0)
            ->u(0.5, 0.6, 0.7, 0);

        $operations = $qc->operations();
        $this->assertSame(['rx', 'ry', 'rz', 'p', 'u'], array_column($operations, 'name'));
        foreach ($operations as $operation) {
            $this->assertInstanceOf(NDArray::class, $operation['params']);
            $this->assertSame([2, 2], $operation['params']->shape());
        }
    }

    public function testEmptyQubitListIsRejected(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new QuantumCircuit(2))->h([]);
    }

    public function testAppendStringGate(): void
    {
        $qc = new QuantumCircuit(2);
        $qc->append('RZ', [1], [0.25]);

        $operations = $qc->operations();
        $this->assertSame('rz', $operations[0]['name']);
        $this->assertSame([1], $operations[0]['qubits']);
        $this->assertSame([0.25], $operations[0]['params']);
        $this->assertSame([], $operations[0]['clbits']);
    }

    public function testAppendGateObject(): void
    {
        $qc = new QuantumCircuit(2);
        $gate = new HGate();
        $qc->append($gate, [1]);

        $operations = $qc->operations();
        $this->assertSame('h', $operations[0]['name']);
        $this->assertSame([1], $operations[0]['qubits']);
        $this->assertInstanceOf(NDArray::class, $operations[0]['params']);
        $this->assertSame([2, 2], $operations[0]['params']->shape());
    }

    public function testAppendGateObjectValidatesQubitCount(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new QuantumCircuit(2))->append(new HGate(), [0, 1]);
    }

    public function testAppendValidatesQubitAndClbitRange(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new QuantumCircuit(1))->append('x', [1]);
    }

    public function testAppendValidatesClbitRange(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new QuantumCircuit(1))->append('measure', [0], [], [0]);
    }

    public function testTwoQubitBroadcasting(): void
    {
        $qc = new QuantumCircuit(3);
        $qc->cx(0, [1, 2]);
        $this->assertSame([[0, 1], [0, 2]], array_column($qc->operations(), 'qubits'));

        $qc = new QuantumCircuit(3);
        $qc->cnot([0, 1], 2);
        $this->assertSame(['cx', 'cx'], array_column($qc->operations(), 'name'));
        $this->assertSame([[0, 2], [1, 2]], array_column($qc->operations(), 'qubits'));

        $qc = new QuantumCircuit(3);
        $qc->cz([0, 1], [2, 0]);
        $this->assertSame([[0, 2], [1, 0]], array_column($qc->operations(), 'qubits'));

        $qc = new QuantumCircuit(3);
        $qc->cp(0.5, 0, [1, 2]);
        $this->assertSame([[0, 1], [0, 2]], array_column($qc->operations(), 'qubits'));

        $qc = new QuantumCircuit(3);
        $qc->swap([0, 1], [2, 0]);
        $this->assertSame([[0, 2], [1, 0]], array_column($qc->operations(), 'qubits'));
    }

    public function testTwoQubitListSizeMismatch(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new QuantumCircuit(3))->cx([0, 1], [2, 0, 1]);
    }

    public function testControlAndTargetMustDiffer(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new QuantumCircuit(2))->cx(0, 0);
    }

    public function testSwapQubitsMustDiffer(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new QuantumCircuit(2))->swap(1, 1);
    }

    public function testUnitaryHelper(): void
    {
        $qc = new QuantumCircuit(2);
        $qc->unitary([[0, 1], [1, 0]], [0, 1]);

        $this->assertSame(['unitary', 'unitary'], array_column($qc->operations(), 'name'));
        $this->assertSame([[0], [1]], array_column($qc->operations(), 'qubits'));

        $qc = new QuantumCircuit(2);
        $qc->unitary([
            [1, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0],
            [0, 1, 0, 0],
        ], [0, 1]);

        $operations = $qc->operations();
        $this->assertCount(1, $operations);
        $this->assertSame([0, 1], $operations[0]['qubits']);
    }

    public function testMeasureMap(): void
    {
        $qc = new QuantumCircuit(2, 2);
        $qc->measure([0, 1], [1, 0]);

        $this->assertSame([0 => 1, 1 => 0], $qc->measureMap());
        $this->assertSame(['measure', 'measure'], array_column($qc->operations(), 'name'));
        $this->assertSame([[1], [0]], array_column($qc->operations(), 'clbits'));
    }

    public function testMeasureSizeMismatch(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new QuantumCircuit(2, 2))->measure([0, 1], [0]);
    }

    public function testMeasureAll(): void
    {
        $qc = new QuantumCircuit(2);
        $qc->measureAll();

        $this->assertSame(2, $qc->numClbits());
        $this->assertSame([0 => 0, 1 => 1], $qc->measureMap());
        $this->assertCount(2, $qc->operations());

        $qc->measureAll();
        $this->assertCount(2, $qc->operations());
    }

    public function testMeasureAllKeepsExistingMeasurements(): void
    {
        $qc = new QuantumCircuit(2, 2);
        $qc->measure(0, 1)->measureAll();

        $this->assertSame([0 => 1, 1 => 1], $qc->measureMap());
        $this->assertCount(2, $qc->operations());
    }

    public function testEffectiveMeasurementInformation(): void
    {
        $qc = new QuantumCircuit(2);
        $this->assertSame([0 => 0, 1 => 1], $qc->effectiveMeasureMap());
        $this->assertSame(2, $qc->effectiveNumClbits());

        $qc = new QuantumCircuit(2, 3);
        $this->assertSame(3, $qc->effectiveNumClbits());

        $qc->measure(0, 2);
        $this->assertSame([0 => 2], $qc->effectiveMeasureMap());
    }

    public function testBarrier(): void
    {
        $qc = new QuantumCircuit(3);
        $qc->barrier();
        $qc->barrier(0, 2);
        $qc->barrier([0, 1], 2);
        $qc->barrier(new ArrayIterator([1, 0]));

        $operations = $qc->operations();
        $this->assertSame([0, 1, 2], $operations[0]['qubits']);
        $this->assertSame([0, 2], $operations[1]['qubits']);
        $this->assertSame([0, 1, 2], $operations[2]['qubits']);
        $this->assertSame([1, 0], $operations[3]['qubits']);
        $this->assertSame('barrier', $operations[0]['name']);
    }

    public function testBarrierValidatesQubitRange(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new QuantumCircuit(2))->barrier(0, 2);
    }

    public function testCountOps(): void
    {
        $qc = new QuantumCircuit(2);
        $qc->h(0)
            ->cx(0, 1)
            ->barrier()
            ->measureAll();

        $this->assertSame(['h' => 1, 'cx' => 1, 'measure' => 2], $qc->countOps());
    }

    public function testDepth(): void
    {
        $qc = new QuantumCircuit(2);
        $qc->h(0)->cx(0, 1);
        $this->assertSame(2, $qc->depth());

        $qc = new QuantumCircuit(2);
        $qc->h(0)->x(1);
        $this->assertSame(1, $qc->depth());

        $qc = new QuantumCircuit(2);
        $qc->h(0)->x(0);
        $this->assertSame(2, $qc->depth());

        $qc->measureAll()->barrier();
        $this->assertSame(2, $qc->depth());
    }

    public function testDraw(): void
    {
        $qc = new QuantumCircuit(2);
        $qc->h(0)->cx(0, 1)->measureAll();
        $drawing = $qc->draw();

        $this->assertStringContainsString('q0:', $drawing);
        $this->assertStringContainsString('q1:', $drawing);
        $this->assertStringContainsString('h', $drawing);
        $this->assertStringContainsString('●', $drawing);
        $this->assertStringContainsString('⊕', $drawing);
        $this->assertStringContainsString('M(c0)', $drawing);
        $this->assertStringContainsString('M(c1)', $drawing);
        $this->assertStringEndsWith("\n", $drawing);
        $this->assertCount(2, explode("\n", rtrim($drawing, "\n")));

        $this->assertSame($drawing, (string)$qc);
    }

    public function testDrawStringGates(): void
    {
        $qc = new QuantumCircuit(2);
        $qc->append('rz', [0], [0.5])
            ->append('qft', [0, 1])
            ->barrier();

        $drawing = $qc->draw();
        $this->assertStringContainsString('rz(0.5)', $drawing);
        $this->assertStringContainsString('QFT', $drawing);
        $this->assertStringContainsString('‖', $drawing);
    }
}
