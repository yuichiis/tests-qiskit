<?php
require __DIR__ . '/../../vendor/autoload.php';

# 位相推定 (QPE: Quantum Phase Estimation)
# T ゲートの位相（45度）を読み取る例題だよ！
# 固有状態 |1> を用意する
# 制御ユニタリを繰り返して 位相を kickback させる
# 逆フーリエ変換 + 測定で位相が分かるよ
# 学びポイント：
# ショア法や化学計算の心臓部だよ
# 量子ビットを増やすと 精度が上がるよ

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Visualization\Visualization;

// 位相推定の回路を作成
$n = 3;  // 位相を推定する精度 (2^n で割った値が誤差)
$qc = new QuantumCircuit($n + 1, $n);

// 固有状態 |1> を用意 (T|1> = e^{iπ/4}|1>)
$qc->x($n);

// 補助ビットを重ね合わせ
$qc->h(range(0, $n - 1));
$qc->barrier();

// 制御ユニタリ (qubit j が U^{2^j} を制御、U=T)
// T^{2^j} の固有値は e^{iπ2^j/4} なので制御位相でかける
// ここでは T ゲートの位相 φ=π/4 を推定する
for ($j = 0; $j < $n; $j++) {
    $qc->cp(M_PI * 2 ** $j / 4, $j, $n);
}
$qc->barrier();

// 逆量子フーリエ変換 (4-1 の QFT と逆向き・逆符号で、H は下位から)
for ($j = 0; $j < intdiv($n, 2); $j++) {
    $qc->swap($j, $n - 1 - $j);
}
for ($j = 0; $j < $n; $j++) {
    for ($k = $j - 1; $k >= 0; $k--) {
        $qc->cp(-M_PI / 2 ** ($j - $k), $k, $j);
    }
    $qc->h($j);
}
$qc->barrier();

// 測定
$qc->measure(range(0, $n - 1), range(0, $n - 1));
echo $qc->draw();

$sim = new StatevectorSampler();
$result = $sim->run([$qc], shots: 1024)->result();
$counts = $result[0]->getCounts();
Visualization::showHistogram($counts);
# 期待：{'001': 1024}（Tゲートの位相 π/4 = 2π/8 を3ビットで読み取る）
