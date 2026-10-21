<?php
require __DIR__ . '/../../vendor/autoload.php';

# 量子テレポーテーション
# 量子の状態を 別の場所に移す仕組みだよ。通信 + もつれを使うよ！

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Visualization\Visualization;

// 3量子ビット：送りたい状態(0) + もつれペア(1,2)
$qc = new QuantumCircuit(3, 3);

# 送りたい状態を作る（ここでは |1> にしてみる）
$qc->x(0);

# もつれペアを作る
$qc->h(1);
$qc->cx(1, 2);

# 送り側の操作
$qc->cx(0, 1);
$qc->h(0);
$qc->measure(0, 0);
$qc->measure(1, 1);

# 受け側の修正（条件付き）
$qc->cx(1, 2);
$qc->cz(0, 2);
$qc->measure(2, 2);

echo $qc->draw();

$sim = new StatevectorSampler();
$result = $sim->run([$qc], shots: 512)->result();
$counts = $result[0]->getCounts();
Visualization::showHistogram($counts);
# 期待：受け側(右端)がほぼいつも `1` になる
