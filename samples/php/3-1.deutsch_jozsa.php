<?php
require __DIR__ . '/../../vendor/autoload.php';

# ドイチュ・ジョサのアルゴリズム
# 関数の性質を 1回の質問で見分ける例題だよ

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Visualization\Visualization;

# 平衡関数の例（CXで作る）
$qc = new QuantumCircuit(2, 1);
$qc->h([0, 1]);
$qc->cx(0, 1);
$qc->h(0);
$qc->measure(0, 0);
print($qc->draw());

$sim = new StatevectorSampler();
$result = $sim->run([$qc], shots: 512)->result();
$counts = $result[0]->getCounts();
Visualization::showHistogram($counts);
# 期待：平衡関数なら `1` が出る、定数関数なら `0` が出る
