<?php
require __DIR__ . '/../../vendor/autoload.php';
# グローバー探索（4件から1件探し）

# 4つの中から当たりを効率よく探す例題だよ。2量子ビットで手計算できるのが嬉しいね

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Visualization\Visualization;

# 当たりを |11> にする例
$qc = new QuantumCircuit(2, 2);
$qc->h([0, 1]);

# オラクル：|11> に印をつける
$qc->cz(0, 1);

# 拡散操作
$qc->h([0, 1]);
$qc->x([0, 1]);
$qc->cz(0, 1);
$qc->x([0, 1]);
$qc->h([0, 1]);

$qc->measure([0, 1], [0, 1]);
echo $qc->draw();

$sim = new StatevectorSampler();
$result = $sim->run([$qc], shots: 512)->result();
$counts = $result[0]->data->meas->getCounts();
Visualization::showHistogram($counts);
# 期待：{'11': ほぼ500} になる
