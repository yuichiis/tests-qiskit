<?php
require __DIR__ . '/../../vendor/autoload.php';

# ベル状態
# 2量子ビットのもつれだよ。測ると 00 か 11 しか出ないのが不思議だね

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Visualization\Visualization;

$qc = new QuantumCircuit(2, 2);
$qc->h(0);
$qc->cx(0, 1);
$qc->measure([0, 1], [0, 1]);
echo $qc->draw();

$sim = new StatevectorSampler();
$result = $sim->run([$qc], shots: 1024)->result();
$counts = $result[0]->getCounts();
Visualization::showHistogram($counts);
# 期待：{'00': 約500, '11': 約500}、`01` と `10` はほぼ0
