<?php
require __DIR__ . '/../../vendor/autoload.php';

# グローバー探索（4件から1件探し）MCXGate 版

# 4つの中から当たりを効率よく探す例題だよ。2量子ビットで手計算できるのが嬉しいね
# MCXGate を使うと、複数の制御ビットを持つ X ゲートをまとめて書けるよ！

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Circuit\MCXGate;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Visualization\Visualization;

# 複数制御の Z ゲート（位相反転）を MCXGate で作る
function mcz(QuantumCircuit $qc, array $qubits): void
{
    $target = $qubits[count($qubits) - 1];
    $qc->h($target);
    $qc->append(new MCXGate(count($qubits) - 1), $qubits);
    $qc->h($target);
}

# 当たりを |11> にする例
$qc = new QuantumCircuit(2, 2);
$qc->h([0, 1]);

# オラクル：|11> に印をつける
mcz($qc, [0, 1]);

# 拡散操作
$qc->h([0, 1]);
$qc->x([0, 1]);
mcz($qc, [0, 1]);
$qc->x([0, 1]);
$qc->h([0, 1]);

$qc->measure([0, 1], [0, 1]);
echo $qc->draw();

$sim = new StatevectorSampler();
$result = $sim->run([$qc], shots: 500)->result();
$counts = $result[0]->getCounts();
Visualization::showHistogram($counts);
# 期待：{'11': ほぼ500} になる
