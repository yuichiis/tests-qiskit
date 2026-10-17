<?php
require __DIR__ . '/../../vendor/autoload.php';

# ショア法（15 = 3 × 5 の因数分解・概念編）
# 本物のショア法は重いので シミュレーションでは考え方だけ掴もうね！
# 15の因数分解は a=7 を選ぶと 7^4 ≡ 1 mod 15 の周期4が見つかるよ
# 周期が分かれば gcd(a^{r/2} ± 1, N) で因数が割れるよ
# 学びポイント：
# 周期探し → 因数分解の流れが大事だよ
# 実機・istsim でも 15や21の小さい数から試そうね

use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Circuit\UnitaryGate;
use Rindow\Qubit\Primitives\StatevectorSampler;
use Rindow\Qubit\Visualization\Visualization;

function shor_gcd(int $a, int $b): int
{
    $a = abs($a);
    $b = abs($b);
    while ($b != 0) {
        [$a, $b] = [$b, $a % $b];
    }
    return $a;
}

function shor_powmod(int $base, int $exp, int $mod): int
{
    $r = 1 % $mod;
    $base %= $mod;
    while ($exp > 0) {
        if ($exp & 1) {
            $r = ($r * $base) % $mod;
        }
        $base = ($base * $base) % $mod;
        $exp >>= 1;
    }
    return $r;
}

/** 連分数で分母上限つきの近似分数 [p, q] を返す */
function shor_limit_denominator(int $num, int $den, int $maxDen): array
{
    $cf = [];
    while ($den != 0) {
        $cf[] = intdiv($num, $den);
        [$num, $den] = [$den, $num % $den];
    }
    $p0 = 1;
    $p1 = $cf[0];
    $q0 = 0;
    $q1 = 1;
    $best = [$p1, $q1];
    for ($i = 1; $i < count($cf); $i++) {
        $p2 = $cf[$i] * $p1 + $p0;
        $q2 = $cf[$i] * $q1 + $q0;
        if ($q2 > $maxDen) {
            break;
        }
        $best = [$p2, $q2];
        [$p0, $p1] = [$p1, $p2];
        [$q0, $q1] = [$q1, $q2];
    }
    return $best;
}

$N = 15;
$a = 7;
$n = 3;  // カウント用量子ビット数 (周期4は 3ビットでちょうど表せる)
$m = 4;  // N を表すのに必要なビット数 (2^4 = 16 > 15)

$qc = new QuantumCircuit($n + $m, $n);

// 作業レジスタを |1> に (モジュラ指数関数の始点)
$qc->x($n);

// カウント側を重ね合わせ
$qc->h(range(0, $n - 1));
$qc->barrier();

// 制御モジュラ乗算 U^{2^j} (|y> -> |a^{2^j} y mod N>)
// j=0: ×7、j=1: ×7^2 ≡ ×4、j=2: ×7^4 ≡ ×1 (恒等変換なので省略)
// ×4 mod 15 は qubitペアの入れ替えと等価なので Fredkin で書ける
$qc->append(UnitaryGate::fredkin(), [1, $n + 0, $n + 2]);
$qc->append(UnitaryGate::fredkin(), [1, $n + 1, $n + 3]);

// ×7 mod 15 は置換ユニタリで作る (|15> は使わないので固定)
$perm = array_fill(0, 16, array_fill(0, 16, 0.0));
for ($x = 0; $x < 15; $x++) {
    $perm[($a * $x) % $N][$x] = 1.0;
}
$perm[15][15] = 1.0;
$qc->append((new UnitaryGate($perm))->control(), [0, $n + 0, $n + 1, $n + 2, $n + 3]);
$qc->barrier();

// 逆量子フーリエ変換 (カウント側、4-1 の QFT と逆向き・逆符号で)
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
$counts = $result[0]->data->meas->getCounts();
Visualization::showHistogram($counts);

// 古典後処理：連分数で周期 r を求めて因数にする
$found = [];
foreach ($counts as $bits => $c) {
    $s = str_pad((string)$bits, $n, '0', STR_PAD_LEFT);
    [, $r] = shor_limit_denominator(bindec($s), 2 ** $n, $N);
    for ($t = $r; $t < $N; $t += $r) {  // r が真の周期の約数だった場合も倍数で試す
        if ($t % 2 == 0 && shor_powmod($a, $t, $N) == 1) {
            $e = $a ** intdiv($t, 2);
            foreach ([shor_gcd($e - 1, $N), shor_gcd($e + 1, $N)] as $f) {
                if ($f > 1 && $f < $N) {
                    $found[$f] = true;
                }
            }
            break;
        }
    }
}
ksort($found);
echo '見つかった因数: [' . implode(', ', array_keys($found)) . "]\n";
# 期待：測定は 000/010/100/110 の4通り、因数は [3, 5]
