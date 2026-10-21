<?php
require __DIR__ . '/../../vendor/autoload.php';

# 2量子ビット位相推定 (QPE: Quantum Phase Estimation)

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Circuit\MCPhaseGate;
use Rindow\Qubit\Circuit\QFTGate;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Visualization\Visualization;

$n = 3;                # 位相を読み取る精度ビット数 (qubit 0..2)
$t0 = $n;              # 対象の2量子ビット (qubit 3, 4)
$t1 = $n + 1;
$phi = M_PI / 4;       # U は |11> にだけ e^{iφ} を掛ける

$qc = new QuantumCircuit($n + 2, $n);

# 固有状態 |11> を用意 (U|11> = e^{iφ}|11>)
$qc->x([$t0, $t1]);

# 補助ビットを重ね合わせ
$qc->h(range(0, $n - 1));
$qc->barrier();

# 制御ユニタリ: qubit j が U^{2^j} を制御
# 制御2個 (補助ビット j と t0)、ターゲット t1 の MCPhaseGate
for ($j = 0; $j < $n; $j++) {
    $qc->append(new MCPhaseGate($phi * 2 ** $j, 2), [$j, $t0, $t1]);
}
$qc->barrier();

# 逆量子フーリエ変換 (QFTGate を使用)
$qc->append((new QFTGate($n))->inverse(), range(0, $n - 1));
$qc->barrier();

$qc->measure(range(0, $n - 1), range(0, $n - 1));
echo $qc->draw();

$sim = new StatevectorSampler();
$result = $sim->run([$qc], shots: 1024)->result();
$counts = $result[0]->getCounts();
Visualization::showHistogram($counts);
# 期待：φ/2π = 1/8 なので 3ビットで '001' がほぼ 1024
