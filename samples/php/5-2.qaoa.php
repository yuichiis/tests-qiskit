<?php
require __DIR__ . '/../../vendor/autoload.php';

# 量子計算 QAOA（Max-Cut・4点リング）

# グラフを2つに分けて 切れる辺を最大化する例題だよ。組み合わせ最適化の入口だよ！
# QAOAは量子回路のパラメータ(beta, gamma)を古典オプティマイザで
# 調整してカット数の期待値を最大化するハイブリッド算法だよ！

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Quantum\StatevectorSimulator;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Visualization\Visualization;
use Rindow\Qubit\Math\NelderMead;

// 4点リング：辺 (0-1, 1-2, 2-3, 3-0)
$edges = [[0, 1], [1, 2], [2, 3], [3, 0]];
$p = 1;

function build_qaoa(float $gamma, float $beta, array $edges): QuantumCircuit
{
    $qc = new QuantumCircuit(4);
    $qc->h(range(0, 3));

    // 問題ユニタリ（ZZ相互作用）
    foreach ($edges as [$i, $j]) {
        $qc->cx($i, $j);
        $qc->rz(2 * $gamma, $j);
        $qc->cx($i, $j);
    }

    // ミキサー
    $qc->rx(2 * $beta, range(0, 3));
    return $qc;
}

function cut_value(string $bits, array $edges): int
{
    // ビット列は左が上位qubitなので右端が q0
    $b = [];
    for ($i = 0; $i < 4; $i++) {
        $b[$i] = (int)$bits[strlen($bits) - 1 - $i];
    }
    $cuts = 0;
    foreach ($edges as [$i, $j]) {
        if ($b[$i] != $b[$j]) {
            $cuts++;
        }
    }
    return $cuts;
}

function neg_expected_cut(array $params, array $edges): float
{
    // 状態ベクトルからカット数の期待値を正確に求める (最適化用)
    [$gamma, $beta] = $params;
    $qc = build_qaoa($gamma, $beta, $edges);
    $sim = new StatevectorSimulator($qc->numQubits());
    $sim->applyCircuit($qc);
    $exp = 0.0;
    foreach ($sim->probabilities() as $bits => $prob) {
        $exp += $prob * cut_value($bits, $edges);
    }
    return -$exp;
}

// 古典オプティマイザでパラメータを最適化 (期待カット数最大化)
[$opt, $fopt] = NelderMead::minimize(
    static fn(array $params) => neg_expected_cut($params, $edges),
    [0.5, 0.5]
);
[$gamma, $beta] = $opt;
printf("最適パラメータ: gamma=%.4f, beta=%.4f (期待カット数 %.4f)\n", $gamma, $beta, -$fopt);

// 最適パラメータで測定
$qc = build_qaoa($gamma, $beta, $edges);
$qc->measureAll();
echo $qc->draw();

$sim = new StatevectorSampler();
$result = $sim->run([$qc], shots: 1000)->result();
$counts = $result[0]->data->meas->getCounts();
Visualization::showHistogram($counts);

// 答えの表示：カット数（両端が別グループの辺の数）が最大のビット列を探す
// (カウントのキーは数値化されることがあるので4桁に正規化する)
$norm = static fn($bits) => str_pad((string)$bits, 4, '0', STR_PAD_LEFT);
$best = 0;
foreach ($counts as $bits => $c) {
    $best = max($best, cut_value($norm($bits), $edges));
}
$winners = [];
foreach ($counts as $bits => $c) {
    if (cut_value($norm($bits), $edges) == $best) {
        $winners[] = $norm($bits);
    }
}
sort($winners);
echo '答え: [' . implode(', ', $winners) . "] (カット数 {$best} / 全4辺)\n";
foreach ($winners as $bits) {
    $c = $counts[$bits] ?? 0;
    echo "  {$bits}: {$c}回\n";
}
