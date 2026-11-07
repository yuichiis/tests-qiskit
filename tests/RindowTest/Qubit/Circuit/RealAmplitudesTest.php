<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Circuit;

use InvalidArgumentException;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Circuit\RealAmplitudes;
use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\QuantumInfo\Statevector;

class RealAmplitudesTest extends TestCase
{
    public function testParameterCountAndStructure(): void
    {
        $ansatz = new RealAmplitudes(numQubits: 2, reps: 2, entanglement: 'linear');

        $this->assertSame(2, $ansatz->numQubits());
        $this->assertSame(2, $ansatz->num_qubits());
        $this->assertSame(2, $ansatz->reps());
        $this->assertSame(6, $ansatz->numParameters());
        $this->assertSame(6, $ansatz->num_parameters());

        $operations = $ansatz->decompose()->operations();
        $this->assertSame(
            ['ry', 'ry', 'cx', 'ry', 'ry', 'cx', 'ry', 'ry'],
            array_column($operations, 'name')
        );
        $this->assertSame(
            [[0], [1], [0, 1], [0], [1], [0, 1], [0], [1]],
            array_column($operations, 'qubits')
        );
    }

    public function testParametrizedDraw(): void
    {
        $ansatz = new RealAmplitudes(numQubits: 2, reps: 1, entanglement: 'linear');
        $drawing = $ansatz->draw();

        $this->assertStringContainsString('ry(θ[0])', $drawing);
        $this->assertStringContainsString('ry(θ[3])', $drawing);
        $this->assertStringContainsString('●', $drawing);
        $this->assertSame($drawing, (string)$ansatz);
    }

    public function testAssignParameters(): void
    {
        $ansatz = new RealAmplitudes(numQubits: 2, reps: 2, entanglement: 'linear');
        $theta = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        $circuit = $ansatz->assignParameters($theta);
        $bind = $ansatz->bind($theta);

        $this->assertSame($ansatz->numQubits(), $circuit->numQubits());
        $this->assertSame(
            ['ry', 'ry', 'cx', 'ry', 'ry', 'cx', 'ry', 'ry'],
            array_column($circuit->operations(), 'name')
        );
        $this->assertSame(
            [[0], [1], [0, 1], [0], [1], [0, 1], [0], [1]],
            array_column($circuit->operations(), 'qubits')
        );

        // 手組みした回路と同じゲート列・行列になることを確認する。
        $expected = (new QuantumCircuit(2))
            ->ry(0.1, 0)
            ->ry(0.2, 1)
            ->cx(0, 1)
            ->ry(0.3, 0)
            ->ry(0.4, 1)
            ->cx(0, 1)
            ->ry(0.5, 0)
            ->ry(0.6, 1);
        $this->assertEquals($expected->operations(), $circuit->operations());
        $this->assertEquals($circuit->operations(), $bind->operations());
    }

    public function testWrongParameterCount(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new RealAmplitudes(2, 2, 'linear'))->assignParameters([0.1, 0.2]);
    }

    public function testZeroReps(): void
    {
        $ansatz = new RealAmplitudes(numQubits: 3, reps: 0);

        $this->assertSame(3, $ansatz->numParameters());
        $this->assertSame(['ry', 'ry', 'ry'], array_column($ansatz->decompose()->operations(), 'name'));
    }

    public function testEntanglementOptions(): void
    {
        $full = new RealAmplitudes(3, 1, 'full');
        $linear = new RealAmplitudes(3, 1, 'linear');
        $reverse = new RealAmplitudes(3, 1, 'reverse_linear');
        $pairwise = new RealAmplitudes(3, 1, 'pairwise');
        $custom = new RealAmplitudes(3, 1, [[0, 2]]);

        $this->assertSame([[0, 1], [0, 2], [1, 2]], $full->entanglementPairs());
        $this->assertSame([[0, 1], [1, 2]], $linear->entanglementPairs());
        $this->assertSame([[1, 0], [2, 1]], $reverse->entanglementPairs());
        $this->assertSame([[0, 1], [1, 2], [2, 0]], $pairwise->entanglementPairs());
        $this->assertSame([[0, 2]], $custom->entanglementPairs());
    }

    public function testInvalidEntanglement(): void
    {
        $this->expectException(InvalidArgumentException::class);
        new RealAmplitudes(2, 1, 'bogus');
    }

    public function testAllZeroParametersGiveZeroState(): void
    {
        $ansatz = new RealAmplitudes(2, 2, 'linear');
        $circuit = $ansatz->assignParameters(array_fill(0, 6, 0.0));
        $probabilities = (new Statevector($circuit))->probabilities();

        $this->assertSame([1.0, 0.0, 0.0, 0.0], $probabilities);
    }
}
