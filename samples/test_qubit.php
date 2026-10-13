<?php
// test.py のPHP版 (Rindow\Qubit 使用)
require __DIR__ . '/../vendor/autoload.php';

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Visualization\Visualization;

echo "hello\n";

// 2. Build your circuit (Bell state)
$qc = new QuantumCircuit(2);
$qc->h(0);
$qc->cx(0, 1);
$qc->measureAll();

echo $qc->draw();

$sampler = new StatevectorSampler();
$result = $sampler->run([$qc], shots: 1024)->result();
$counts = $result[0]->data->meas->getCounts();
//print_r($counts);
// 期待値: '00' と '11' がほぼ半々

// 3. Visualize your results
Visualization::showHistogram($counts);

// 4. See what happens: 3qubit GHZ
$qc3 = new QuantumCircuit(3);
$qc3->h(0);
$qc3->cx(0, 1);
$qc3->cx(1, 2);
$qc3->measureAll();
$counts3 = $sampler->run([$qc3], shots: 1024)->result()[0]->data->meas->getCounts();
//print_r($counts3);
Visualization::showHistogram($counts3);
