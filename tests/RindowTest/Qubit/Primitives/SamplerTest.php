<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Primitives;

use InvalidArgumentException;
use LogicException;
use OutOfBoundsException;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Primitives\MeasData;
use Rindow\Qubit\Primitives\PubResult;
use Rindow\Qubit\Primitives\SamplerJob;
use Rindow\Qubit\Primitives\SamplerResult;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\QuantumCircuit;

class SamplerTest extends TestCase
{
    public function testMeasData(): void
    {
        $meas = new MeasData(['00' => 3, '11' => 1], 4);

        $this->assertSame(['00' => 3, '11' => 1], $meas->getCounts());
        $this->assertSame(4, $meas->numShots());
        $this->assertEqualsWithDelta(['00' => 0.75, '11' => 0.25], $meas->getProbabilities(), 1.0e-12);
    }

    public function testRunSingleCircuit(): void
    {
        $qc = (new QuantumCircuit(1))->x(0)->measureAll();
        $sampler = new StatevectorSampler(seed: 123);
        $result = $sampler->run([$qc], shots: 50)->result();

        $this->assertInstanceOf(SamplerResult::class, $result);
        $this->assertCount(1, $result);
        $this->assertInstanceOf(PubResult::class, $result[0]);
        $this->assertSame([1 => 50], $result[0]->getCounts());
        $this->assertSame(50, $result[0]->meas->numShots());
        $this->assertSame(2, $result->getStatevector()->dim());
        $this->assertSame($result[0]->getStatevector(), $result->getStatevector(0));
    }

    public function testRunMultipleCircuits(): void
    {
        $x = (new QuantumCircuit(1))->x(0)->measureAll();
        $bell = (new QuantumCircuit(2))->h(0)->cx(0, 1)->measureAll();
        $sampler = new StatevectorSampler(seed: 1);
        $result = $sampler->run([$x, $bell], shots: 100)->result();

        $this->assertCount(2, $result);
        $this->assertSame([1 => 100], $result[0]->getCounts());

        $bellCounts = $result[1]->getCounts();
        $this->assertSame(100, array_sum($bellCounts));
        foreach (array_keys($bellCounts) as $bits) {
            $this->assertContains((string)$bits, ['00', '11']);
        }
    }

    public function testSeedMakesResultsReproducible(): void
    {
        $qc = (new QuantumCircuit(2))->h(0)->cx(0, 1)->measureAll();
        $sampler = new StatevectorSampler(seed: 9876);

        $first = $sampler->run([$qc], shots: 100)->result()[0]->getCounts();
        $second = $sampler->run([$qc], shots: 100)->result()[0]->getCounts();

        $this->assertSame($first, $second);
    }

    public function testExplicitSeedOverridesDefault(): void
    {
        $qc = (new QuantumCircuit(2))->h(0)->cx(0, 1)->measureAll();

        $first = (new StatevectorSampler(seed: 1))
            ->run([$qc], shots: 100, seed: 42)
            ->result()[0]->getCounts();
        $second = (new StatevectorSampler(seed: 2))
            ->run([$qc], shots: 100, seed: 42)
            ->result()[0]->getCounts();

        $this->assertSame($first, $second);
    }

    public function testRunWithoutExplicitMeasurement(): void
    {
        $qc = (new QuantumCircuit(1))->x(0);
        $result = (new StatevectorSampler(seed: 3))->run([$qc], shots: 10)->result();

        $this->assertSame([1 => 10], $result[0]->getCounts());
    }

    public function testInvalidShots(): void
    {
        $this->expectException(InvalidArgumentException::class);
        new SamplerJob([new QuantumCircuit(1)], 0);
    }

    public function testEmptyCircuits(): void
    {
        $this->expectException(InvalidArgumentException::class);
        new SamplerJob([]);
    }

    public function testInvalidCircuitElement(): void
    {
        $this->expectException(InvalidArgumentException::class);
        new SamplerJob([new \stdClass()]);
    }

    public function testSamplerResultArrayAccessAndIteration(): void
    {
        $qc = (new QuantumCircuit(1))->x(0)->measureAll();
        $result = (new StatevectorSampler(seed: 1))->run([$qc, $qc], shots: 10)->result();

        $this->assertTrue(isset($result[0]));
        $this->assertTrue($result->offsetExists(1));
        $this->assertFalse($result->offsetExists(2));
        $this->assertSame($result->get(0), $result[0]);
        $this->assertSame($result->get(1), $result->offsetGet(1));

        $pubs = [];
        foreach ($result as $pub) {
            $pubs[] = $pub;
        }
        $this->assertCount(2, $pubs);
        $this->assertSame($result->get(0), $pubs[0]);
    }

    public function testSamplerResultIsImmutable(): void
    {
        $result = (new StatevectorSampler())->run([new QuantumCircuit(1)], shots: 1)->result();

        $this->expectException(LogicException::class);
        $result[0] = $result[0];
    }

    public function testSamplerResultCannotUnset(): void
    {
        $result = (new StatevectorSampler())->run([new QuantumCircuit(1)], shots: 1)->result();

        $this->expectException(LogicException::class);
        unset($result[0]);
    }

    public function testSamplerResultOutOfRange(): void
    {
        $result = (new StatevectorSampler())->run([new QuantumCircuit(1)], shots: 1)->result();

        $this->expectException(OutOfBoundsException::class);
        $result->get(1);
    }
}
