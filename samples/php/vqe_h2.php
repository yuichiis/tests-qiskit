<?php
require __DIR__ . '/../../vendor/autoload.php';

/**
 * VQE による H2 分子の基底状態エネルギー計算 (PHP版)
 *
 * - ハミルトニアン: H2, STO-3G, 結合長 0.735 Å を Jordan-Wigner + 2量子ビット縮約したもの
 *   (PySCF なしで動くよう係数を直書き)
 * - ansatz: RealAmplitudes (RY + CX)
 * - 期待値計算: StatevectorEstimator (厳密な状態ベクトル、ショットノイズなし)
 * - オプティマイザ: COBYLA
 * - 精度: Precision の既定値 (float64 + complex128)。
 *   float32 + complex64 で実行するには、計算前に
 *   \Rindow\Qubit\Math\Precision::setFloat32() を呼ぶ。
 */

use Rindow\Qubit\Circuit\RealAmplitudes;
use Rindow\Qubit\Math\HermitianEigen;
use Rindow\Qubit\Math\Minimize;
use Rindow\Qubit\Primitives\StatevectorEstimator;
use Rindow\Qubit\QuantumInfo\SparsePauliOp;

// ---------------------------------------------------------------------------
// 1. ハミルトニアン (電子部分のみ)。核間反発エネルギーは別に足す。
// ---------------------------------------------------------------------------
$nuclearRepulsion = 0.7199689944489797;  // Hartree (0.735 Å)

$hamiltonian = SparsePauliOp::fromList([
    ['II', -1.052373245772859],
    ['IZ', 0.39793742484318045],
    ['ZI', -0.39793742484318045],
    ['ZZ', -0.01128010425623538],
    ['XX', 0.18093119978423156],
]);

// ---------------------------------------------------------------------------
// 2. ansatz (パラメータ付き回路)
// ---------------------------------------------------------------------------
$ansatz = new RealAmplitudes(numQubits: 2, reps: 2, entanglement: 'linear');
echo $ansatz->decompose()->draw();
echo "パラメータ数: {$ansatz->numParameters()}\n\n";

// ---------------------------------------------------------------------------
// 3. コスト関数: <ψ(θ)|H|ψ(θ)>
// ---------------------------------------------------------------------------
$estimator = new StatevectorEstimator();
$history = [];

$cost = static function (array $theta) use ($ansatz, $hamiltonian, $estimator, &$history): float {
    // Estimator の入力は [回路, 観測量, パラメータ値] の PUB
    $pub = [$ansatz, $hamiltonian, $theta];
    $energy = (float)$estimator->run([$pub])->result()[0]->data->evs;
    $history[] = $energy;
    return $energy;
};

// ---------------------------------------------------------------------------
// 4. 古典オプティマイザで最小化
// ---------------------------------------------------------------------------
mt_srand(42);
$theta0 = [];
for ($i = 0; $i < $ansatz->numParameters(); $i++) {
    $theta0[] = mt_rand() / mt_getrandmax() * 2.0 * M_PI - M_PI;
}

$result = Minimize::minimize(
    $cost,
    $theta0,
    method: 'COBYLA',
    options: ['maxiter' => 500, 'tol' => 1.0e-8],
);

$vqeElectronic = $result->fun;
$vqeTotal = $vqeElectronic + $nuclearRepulsion;

// ---------------------------------------------------------------------------
// 5. 厳密対角化との比較 (2量子ビットなので行列を直接作れる)
// ---------------------------------------------------------------------------
// 実行パスも参照値も同じ既定精度 (float64 + complex128) を使う。
$exactElectronic = min(HermitianEigen::eigvalsh(
    $hamiltonian->toMatrix()
));
$exactTotal = $exactElectronic + $nuclearRepulsion;

echo '最適化の反復回数 (関数評価回数): ' . count($history) . "\n";
printf("VQE   電子エネルギー : %.8f Ha\n", $vqeElectronic);
printf("厳密  電子エネルギー : %.8f Ha\n", $exactElectronic);
printf("VQE   全エネルギー   : %.8f Ha\n", $vqeTotal);
printf("厳密  全エネルギー   : %.8f Ha\n", $exactTotal);
printf("誤差                 : %.2e Ha\n", abs($vqeTotal - $exactTotal));
