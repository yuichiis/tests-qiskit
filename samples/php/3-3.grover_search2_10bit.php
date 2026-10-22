<?php
require __DIR__ . '/../../vendor/autoload.php';

# グローバー探索（1024件から1件探し）MCXGate 版

# 10量子ビット（1024通り）の中から当たりを効率よく探す例題だよ
# 当たりが |1111111111> になるようにオラクルを組むよ
# 最適な繰り返し回数は floor(pi/4 * sqrt(1024)) = 25 回！

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Circuit\MCXGate;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Visualization\Visualization;

$N = 10;
$TARGET = str_repeat('1', $N);
$ITERATIONS = (int)floor(pi() / 4 * sqrt(2 ** $N));

# 複数制御の Z ゲート（位相反転）を MCXGate で作る
function mcz(QuantumCircuit $qc, array $qubits): void
{
    $target = $qubits[count($qubits) - 1];
    $qc->h($target);
    $qc->append(new MCXGate(count($qubits) - 1), $qubits);
    $qc->h($target);
}

# 当たりを |111...1> にする例
$qc = new QuantumCircuit($N, $N);
$qc->h(range(0, $N - 1));

# オラクルと拡散操作を最適な回数だけ繰り返す
for ($i = 0; $i < $ITERATIONS; $i++) {
    # オラクル：|111...1> に印をつける
    mcz($qc, range(0, $N - 1));

    # 拡散操作
    $qc->h(range(0, $N - 1));
    $qc->x(range(0, $N - 1));
    mcz($qc, range(0, $N - 1));
    $qc->x(range(0, $N - 1));
    $qc->h(range(0, $N - 1));
}

$qc->measure(range(0, $N - 1), range(0, $N - 1));
echo $qc->draw();
echo "繰り返し回数: {$ITERATIONS}\n";

$sim = new StatevectorSampler();
$result = $sim->run([$qc], shots: 500)->result();
$counts = $result[0]->getCounts();
Visualization::showHistogram($counts);
# 期待：{'1111111111': ほぼ500} になる
