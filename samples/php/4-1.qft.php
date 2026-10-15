<?php
require __DIR__ . '/../../vendor/autoload.php';

# 量子フーリエ変換（2量子ビット）QFT: Quantum Fourier Transform

# 位相推定やショア法の部品だよ。まずは2量子ビットで回してみよう！

use Rindow\Qubit\QuantumCircuit;

$qc = new QuantumCircuit(2);
$qc->h(1);
$qc->cp(M_PI / 2, 0, 1);
$qc->h(0);
$qc->swap(0, 1);
echo $qc->draw();

$state = $qc->statevector();
foreach ($state as $i => $amp) {
    printf("%02b: %.4f%+.4fi\n", $i, $amp['re'], $amp['im']);
}
# 期待：|00>入力のQFTは一様な重ね合わせ (0.5ずつ)
