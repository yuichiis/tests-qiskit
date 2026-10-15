<?php
require __DIR__ . '/../../vendor/autoload.php';

# 量子版コイン投げ
# 2量子ビットのもつれだよ。測ると 00 か 11 しか出ないのが不思議だね

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Visualization\Visualization;


$qc = new QuantumCircuit(1, 1);
$qc->h(0);
$qc->measure(0, 0);
echo $qc->draw();

$sim = new StatevectorSampler();
$result = $sim->run([$qc], shots: 1024)->result();
$counts = $result[0]->data->meas->getCounts();
Visualization::showHistogram($counts);
# 期待：{'0': 約500, '1': 約500}
