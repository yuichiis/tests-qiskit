<?php
require __DIR__ . '/../../vendor/autoload.php';

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Visualization\Visualization;

# バーンスタイン・ヴァジラーニのアルゴリズム
# 隠れたビット列 s を 1回で見つける例題だよ！

# 隠れた文字列 s = 11 の例
$s = "11";
$n = strlen($s);
$qc = new QuantumCircuit($n + 1, $n);

$qc->x($n);
$qc->h(range(0, $n));

# オラクル：sに応じてCXを貼る
foreach (array_reverse(str_split($s)) as $i => $bit) {
    if ($bit == "1") {
        $qc->cx($i, $n);
    }
}

$qc->h(range(0, $n - 1));
$qc->measure(range(0, $n - 1), range(0, $n - 1));
echo $qc->draw();

$sim = new StatevectorSampler();
$result = $sim->run([$qc], shots: 512)->result();
$counts = $result[0]->getCounts();
Visualization::showHistogram($counts);
# 期待：{'11': 500} みたいに隠れた文字列がそのまま出る
