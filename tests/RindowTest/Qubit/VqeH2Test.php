<?php
declare(strict_types=1);

namespace RindowTest\Qubit;

use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Circuit\RealAmplitudes;
use Rindow\Qubit\Math\HermitianEigen;
use Rindow\Qubit\Math\Minimize;
use Rindow\Qubit\Primitives\StatevectorEstimator;
use Rindow\Qubit\QuantumInfo\SparsePauliOp;

/**
 * samples/php/vqe_h2.php と同じ手順を小さく再現するエンドツーエンドテスト。
 */
class VqeH2Test extends TestCase
{
    public function testVqeReachesExactH2GroundState(): void
    {
        $nuclearRepulsion = 0.7199689944489797;

        $hamiltonian = SparsePauliOp::fromList([
            ['II', -1.052373245772859],
            ['IZ', 0.39793742484318045],
            ['ZI', -0.39793742484318045],
            ['ZZ', -0.01128010425623538],
            ['XX', 0.18093119978423156],
        ]);

        $ansatz = new RealAmplitudes(numQubits: 2, reps: 2, entanglement: 'linear');
        $estimator = new StatevectorEstimator();

        $cost = static function (array $theta) use ($ansatz, $hamiltonian, $estimator): float {
            return (float)$estimator
                ->run([[$ansatz, $hamiltonian, $theta]])
                ->result()[0]->data->evs;
        };

        // samples/php/vqe_h2.php の mt_srand(42) で生成される初期値
        $theta0 = [
            -0.788287709314266,
            1.863234522144295,
            2.831921542926562,
            -1.989037890326693,
            1.457660906815001,
            1.757350366497756,
        ];

        $result = Minimize::minimize(
            $cost,
            $theta0,
            method: 'COBYLA',
            options: ['maxiter' => 500, 'tol' => 1.0e-8],
        );

        // 実行パスと同じ既定精度 (float64 + complex128) で参照値を求める。
        $exact = min(HermitianEigen::eigvalsh($hamiltonian->toMatrix()));

        $this->assertTrue($result->success);
        $this->assertGreaterThan(0, $result->nfev);
        $this->assertEqualsWithDelta($exact, $result->fun, 1.0e-8);
        $this->assertEqualsWithDelta($exact + $nuclearRepulsion, $result->fun + $nuclearRepulsion, 1.0e-8);
        $this->assertEqualsWithDelta(-1.8572750302023795, $result->fun, 1.0e-8);
        $this->assertEqualsWithDelta(-1.1373060357533998, $result->fun + $nuclearRepulsion, 1.0e-8);
    }
}
