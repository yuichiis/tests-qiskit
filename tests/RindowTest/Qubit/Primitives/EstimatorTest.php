<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Primitives;

use InvalidArgumentException;
use LogicException;
use OutOfBoundsException;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Circuit\RealAmplitudes;
use Rindow\Qubit\Math\HermitianEigen;
use Rindow\Qubit\Primitives\EstimatorPubResult;
use Rindow\Qubit\Primitives\EstimatorResult;
use Rindow\Qubit\Primitives\StatevectorEstimator;
use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\QuantumInfo\SparsePauliOp;

class EstimatorTest extends TestCase
{
    private function z(): SparsePauliOp
    {
        return SparsePauliOp::fromList([['Z', 1.0]]);
    }

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

    public function testSingleQubitExpectationValues(): void
    {
        $estimator = new StatevectorEstimator();
        $zero = new QuantumCircuit(1);
        $one = (new QuantumCircuit(1))->x(0);

        $result = $estimator->run([[$zero, $this->z()], [$one, $this->z()]])->result();

        $this->assertInstanceOf(EstimatorResult::class, $result);
        $this->assertCount(2, $result);
        $this->assertInstanceOf(EstimatorPubResult::class, $result[0]);
        $this->assertEqualsWithDelta(1.0, $result[0]->data->evs, 1.0e-12);
        $this->assertEqualsWithDelta(-1.0, $result[1]->data->evs, 1.0e-12);
        $this->assertSame(0.0, $result[0]->data->stds);
        $this->assertEqualsWithDelta(1.0, $result->getExpectationValue(0), 1.0e-12);
    }

    public function testAssociativePub(): void
    {
        $result = (new StatevectorEstimator())->run([[
            'circuit' => (new QuantumCircuit(1))->x(0),
            'observable' => $this->z(),
        ]])->result();

        $this->assertEqualsWithDelta(-1.0, $result[0]->getExpectationValue(), 1.0e-12);
    }

    public function testH2GroundStateEnergyAtKnownOptimum(): void
    {
        $ansatz = new RealAmplitudes(2, 2, 'linear');
        $theta = [
            2.7149888827301387,
            0.6978178975908345,
            3.318580283519334,
            1.1460080910787043,
            -2.884078832813352,
            2.7515532489999335,
        ];
        $h = $this->h2();
        // 実行パスと同じ既定精度 (float64 + complex128) で参照値を求める。
        $exact = min(HermitianEigen::eigvalsh($h->toMatrix()));

        $result = (new StatevectorEstimator())
            ->run([[$ansatz, $h, $theta]])
            ->result();

        $this->assertEqualsWithDelta($exact, $result[0]->data->evs, 1.0e-10);
        $this->assertSame(4, $result[0]->getStatevector()->dim());
    }

    public function testElectronCountMismatch(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new StatevectorEstimator())->run([[new QuantumCircuit(2), $this->z()]]);
    }

    public function testInvalidObservable(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new StatevectorEstimator())->run([[new QuantumCircuit(1), new \stdClass()]]);
    }

    public function testEmptyPubs(): void
    {
        $this->expectException(InvalidArgumentException::class);
        (new StatevectorEstimator())->run([]);
    }

    public function testResultArrayAccess(): void
    {
        $result = (new StatevectorEstimator())
            ->run([[new QuantumCircuit(1), $this->z()]])
            ->result();

        $this->assertTrue(isset($result[0]));
        $this->assertFalse($result->offsetExists(1));
        $this->assertSame($result->get(0), $result[0]);

        $pubs = [];
        foreach ($result as $pub) {
            $pubs[] = $pub;
        }
        $this->assertCount(1, $pubs);
    }

    public function testResultIsImmutable(): void
    {
        $result = (new StatevectorEstimator())
            ->run([[new QuantumCircuit(1), $this->z()]])
            ->result();

        $this->expectException(LogicException::class);
        $result[0] = $result[0];
    }

    public function testOutOfRange(): void
    {
        $result = (new StatevectorEstimator())
            ->run([[new QuantumCircuit(1), $this->z()]])
            ->result();

        $this->expectException(OutOfBoundsException::class);
        $result->get(1);
    }
}
