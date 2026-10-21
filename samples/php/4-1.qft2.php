<?php
require __DIR__ . '/../../vendor/autoload.php';

# 量子フーリエ変換（2量子ビット）QFT: Quantum Fourier Transform

# 位相推定やショア法の部品だよ。まずは2量子ビットで回してみよう！
# QFTGate を使うと、手書きの H や制御位相ゲートをまとめて書けるよ！

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Circuit\QFTGate;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Math\IntMath;

$qc = new QuantumCircuit(2);
$qc->append(new QFTGate(2), [0, 1]);
echo $qc->draw();

$sim = new StatevectorSampler();
$result = $sim->run([$qc], shots: 1024)->result();
$state = $result->getStatevector();
foreach ($state as $i => $amp) {
    printf("%s: %.4f%+.4fi\n", IntMath::int2bit($i, 2), $amp->real, $amp->imag);
}
