<?php
require __DIR__ . '/../../vendor/autoload.php';

# 例題2-2：超密度符号（2ビットを1量子ビットで送る）
# もつれがあると 1つの量子ビットで 2ビット分の情報を送れるよ！
#
# 送りたい2ビットで操作を変えるよ：00→何もしない、01→Z、10→X、11→XZ
# 受け側で CX + H + 測定 すると元に戻るよ
# 学びポイント：
#
# もつれは通信の資源になるよ
# テレポートとちょうど逆の関係だよ

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Visualization\Visualization;

// Alice と Bob がもつれを共有した状態からスタート
$qc = new QuantumCircuit(2, 2);
$qc->h(0);
$qc->cx(0, 1);
$qc->barrier();

// Alice が送りたいビット（例：10）
$message = "10";

// 00: do nothing
if ($message == "00") {
    // 何もしない
}
// 01: Z
elseif ($message == "01") {
    $qc->z(0);
}
// 10: X
elseif ($message == "10") {
    $qc->x(0);
}
// 11: XZ
elseif ($message == "11") {
    $qc->x(0);
    $qc->z(0);
}

$qc->barrier();

// Bob 側：CX + H + 測定
$qc->cx(0, 1);
$qc->h(0);
$qc->measure([0, 1], [0, 1]);

echo $qc->draw();

$sim = new StatevectorSampler();
$result = $sim->run([$qc], shots: 1024)->result();
$counts = $result[0]->getCounts();
Visualization::showHistogram($counts);
# 期待：送ったメッセージと同じビット列が100%出る（例：10なら {'10': 1024}）
