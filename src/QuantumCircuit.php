<?php
declare(strict_types=1);

namespace Rindow\Qubit;

use InvalidArgumentException;
use Rindow\Qubit\Quantum\StatevectorSimulator;

/**
 * Qiskit互換の量子回路クラス (PHP版 qiskit::QuantumCircuit)
 *
 * Qiskitのエンディアン規約に従う:
 *  - qubit 0 が最下位ビット (LSB)
 *  - ビット列の左端が最上位の古典ビット (例: 2qubitで '01' は q1=0, q0=1)
 *
 * ```php
 * use Rindow\Qubit\QuantumCircuit;
 * $qc = new QuantumCircuit(2);
 * $qc->h(0);
 * $qc->cx(0, 1);
 * $qc->measureAll();
 * ```
 */
class QuantumCircuit
{
    private int $numQubits;
    private int $numClbits;
    /** @var array<int,array{name:string,qubits:array<int>,params:array<float>,clbits:array<int>}> */
    private array $operations = [];
    /** @var array<int,int> qubit => clbit の測定マップ */
    private array $measureMap = [];

    public function __construct(int $numQubits, int $numClbits = 0)
    {
        if ($numQubits < 1 || $numQubits > 30) {
            throw new InvalidArgumentException('numQubits must be in [1, 30].');
        }
        if ($numClbits < 0) {
            throw new InvalidArgumentException('numClbits must be >= 0.');
        }
        $this->numQubits = $numQubits;
        $this->numClbits = $numClbits;
    }

    public function numQubits(): int
    {
        return $this->numQubits;
    }

    public function numClbits(): int
    {
        return $this->numClbits;
    }

    /** @deprecated use numQubits() */
    public function num_qubits(): int
    {
        return $this->numQubits;
    }

    /** @deprecated use numClbits() */
    public function num_clbits(): int
    {
        return $this->numClbits;
    }

    /** @return array<int,array{name:string,qubits:array<int>,params:array<float>,clbits:array<int>}> */
    public function operations(): array
    {
        return $this->operations;
    }

    /** @return array<int,int> */
    public function measureMap(): array
    {
        return $this->measureMap;
    }

    private function checkQubit(int $q): void
    {
        if ($q < 0 || $q >= $this->numQubits) {
            throw new InvalidArgumentException("Qubit index out of range: {$q}");
        }
    }

    private function checkClbit(int $c): void
    {
        if ($c < 0 || $c >= $this->numClbits) {
            throw new InvalidArgumentException("Clbit index out of range: {$c}");
        }
    }

    private function append(string $name, array $qubits, array $params = [], array $clbits = []): static
    {
        foreach ($qubits as $q) {
            $this->checkQubit($q);
        }
        foreach ($clbits as $c) {
            $this->checkClbit($c);
        }
        $this->operations[] = [
            'name' => strtolower($name),
            'qubits' => array_values($qubits),
            'params' => array_values($params),
            'clbits' => array_values($clbits),
        ];
        return $this;
    }

    // ---- 単一量子ビットゲート ----
    // Qiskit互換: qubit引数は int または intの配列/Traversableを受け付け、
    // 配列の場合は各qubitに順にゲートを適用する.

    /**
     * @param int|array<int>|\Traversable $qubit
     * @return array<int> 正規化されたqubitリスト
     */
    private static function normalizeQubitList(int|array|\Traversable $qubit): array
    {
        if ($qubit instanceof \Traversable) {
            $qubit = iterator_to_array($qubit, false);
        }
        if (is_array($qubit)) {
            $list = array_values($qubit);
            if (count($list) === 0) {
                throw new InvalidArgumentException('Qubit list must not be empty.');
            }
            return $list;
        }
        return [$qubit];
    }

    /**
     * @param int|array<int>|\Traversable $clbit
     * @return array<int>
     */
    private static function normalizeClbitList(int|array|\Traversable $clbit): array
    {
        if ($clbit instanceof \Traversable) {
            $clbit = iterator_to_array($clbit, false);
        }
        if (is_array($clbit)) {
            $list = array_values($clbit);
            if (count($list) === 0) {
                throw new InvalidArgumentException('Clbit list must not be empty.');
            }
            return $list;
        }
        return [$clbit];
    }

    /**
     * 2qubitゲート用の引数展開. Qiskit互換のブロードキャスト:
     * - 両方単一 => 1ペア
     * - 片方が単一、もう一方がリスト => 単一を各要素にブロードキャスト
     * - 両方リストで同数 => 1対1対応
     * - それ以外 => エラー
     * @return array<array{int,int}>
     */
    private static function expandTwoQubitArgs(int|array|\Traversable $q1, int|array|\Traversable $q2): array
    {
        $list1 = self::normalizeQubitList($q1);
        $list2 = self::normalizeQubitList($q2);
        $n1 = count($list1);
        $n2 = count($list2);
        if ($n1 === 1 && $n2 === 1) {
            return [[(int)$list1[0], (int)$list2[0]]];
        }
        if ($n1 === 1) {
            $pairs = [];
            foreach ($list2 as $b) {
                $pairs[] = [(int)$list1[0], (int)$b];
            }
            return $pairs;
        }
        if ($n2 === 1) {
            $pairs = [];
            foreach ($list1 as $a) {
                $pairs[] = [(int)$a, (int)$list2[0]];
            }
            return $pairs;
        }
        if ($n1 !== $n2) {
            throw new InvalidArgumentException(
                "Qubit list size mismatch: {$n1} vs {$n2}."
            );
        }
        $pairs = [];
        foreach ($list1 as $i => $a) {
            $pairs[] = [(int)$a, (int)$list2[$i]];
        }
        return $pairs;
    }

    public function h(int|array|\Traversable $qubit): static
    {
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('h', [(int)$q]);
        }
        return $this;
    }

    public function x(int|array|\Traversable $qubit): static
    {
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('x', [(int)$q]);
        }
        return $this;
    }

    public function y(int|array|\Traversable $qubit): static
    {
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('y', [(int)$q]);
        }
        return $this;
    }

    public function z(int|array|\Traversable $qubit): static
    {
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('z', [(int)$q]);
        }
        return $this;
    }

    public function s(int|array|\Traversable $qubit): static
    {
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('s', [(int)$q]);
        }
        return $this;
    }

    public function sdg(int|array|\Traversable $qubit): static
    {
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('sdg', [(int)$q]);
        }
        return $this;
    }

    public function t(int|array|\Traversable $qubit): static
    {
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('t', [(int)$q]);
        }
        return $this;
    }

    public function tdg(int|array|\Traversable $qubit): static
    {
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('tdg', [(int)$q]);
        }
        return $this;
    }

    public function id(int|array|\Traversable $qubit): static
    {
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('id', [(int)$q]);
        }
        return $this;
    }

    public function rx(float $theta, int|array|\Traversable $qubit): static
    {
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('rx', [(int)$q], [$theta]);
        }
        return $this;
    }

    public function ry(float $theta, int|array|\Traversable $qubit): static
    {
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('ry', [(int)$q], [$theta]);
        }
        return $this;
    }

    public function rz(float $theta, int|array|\Traversable $qubit): static
    {
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('rz', [(int)$q], [$theta]);
        }
        return $this;
    }

    /** Phase gate P(λ) == U1(λ) */
    public function p(float $lambda, int|array|\Traversable $qubit): static
    {
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('p', [(int)$q], [$lambda]);
        }
        return $this;
    }

    public function u(float $theta, float $phi, float $lambda, int|array|\Traversable $qubit): static
    {
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('u', [(int)$q], [$theta, $phi, $lambda]);
        }
        return $this;
    }

    // ---- 2量子ビットゲート ----

    public function cx(int|array|\Traversable $control, int|array|\Traversable $target): static
    {
        foreach (self::expandTwoQubitArgs($control, $target) as [$c, $t]) {
            if ($c === $t) {
                throw new InvalidArgumentException('Control and target must differ.');
            }
            $this->append('cx', [$c, $t]);
        }
        return $this;
    }

    /** alias of cx() */
    public function cnot(int|array|\Traversable $control, int|array|\Traversable $target): static
    {
        return $this->cx($control, $target);
    }

    public function cz(int|array|\Traversable $control, int|array|\Traversable $target): static
    {
        foreach (self::expandTwoQubitArgs($control, $target) as [$c, $t]) {
            if ($c === $t) {
                throw new InvalidArgumentException('Control and target must differ.');
            }
            $this->append('cz', [$c, $t]);
        }
        return $this;
    }

    public function swap(int|array|\Traversable $qubit1, int|array|\Traversable $qubit2): static
    {
        foreach (self::expandTwoQubitArgs($qubit1, $qubit2) as [$a, $b]) {
            if ($a === $b) {
                throw new InvalidArgumentException('Swap qubits must differ.');
            }
            $this->append('swap', [$a, $b]);
        }
        return $this;
    }

    /**
     * 任意の1qubitユニタリ行列 [[a,b],[c,d]] を適用する.
     * 各要素は float(実数) または [re, im] / ['re'=>..,'im'=>..] / Complex オブジェクト.
     */
    public function unitary(array $matrix, int|array|\Traversable $qubit): static
    {
        $flat = self::matrixToComplexPairs($matrix);
        if (count($flat) !== 4) {
            throw new InvalidArgumentException('unitary() requires a 2x2 matrix.');
        }
        foreach (self::normalizeQubitList($qubit) as $q) {
            $this->append('unitary', [(int)$q], $flat);
        }
        return $this;
    }

    /**
     * @param array $matrix 2x2
     * @return array<float> [re00,im00, re01,im01, re10,im10, re11,im11]
     */
    private static function matrixToComplexPairs(array $matrix): array
    {
        $rows = array_values($matrix);
        if (count($rows) !== 2 || count((array)$rows[0]) !== 2 || count((array)$rows[1]) !== 2) {
            throw new InvalidArgumentException('unitary() requires a 2x2 matrix.');
        }
        $out = [];
        foreach ($rows as $row) {
            foreach ((array)$row as $v) {
                [$re, $im] = self::toComplexPair($v);
                $out[] = $re;
                $out[] = $im;
            }
        }
        return $out;
    }

    /** @return array{float,float} */
    private static function toComplexPair(mixed $v): array
    {
        if (is_array($v)) {
            if (array_key_exists('re', $v) || array_key_exists('im', $v)) {
                return [(float)($v['re'] ?? 0.0), (float)($v['im'] ?? 0.0)];
            }
            if (count($v) === 2) {
                $a = array_values($v);
                return [(float)$a[0], (float)$a[1]];
            }
        }
        if (is_object($v)) {
            if (isset($v->real) || isset($v->imag)) {
                return [(float)($v->real ?? 0.0), (float)($v->imag ?? 0.0)];
            }
            if (isset($v->re) || isset($v->im)) {
                return [(float)($v->re ?? 0.0), (float)($v->im ?? 0.0)];
            }
        }
        return [(float)$v, 0.0];
    }

    // ---- 測定・バリア ----

    public function measure(int|array|\Traversable $qubit, int|array|\Traversable $clbit): static
    {
        $qubits = self::normalizeQubitList($qubit);
        $clbits = self::normalizeClbitList($clbit);
        if (count($qubits) !== count($clbits)) {
            throw new InvalidArgumentException(
                'measure() qubit/clbit size mismatch: ' . count($qubits) . ' vs ' . count($clbits) . '.'
            );
        }
        foreach ($qubits as $i => $q) {
            $q = (int)$q;
            $c = (int)$clbits[$i];
            $this->checkQubit($q);
            $this->checkClbit($c);
            $this->measureMap[$q] = $c;
            $this->append('measure', [$q], [], [$c]);
        }
        return $this;
    }

    public function measureAll(): static
    {
        return $this->measure_all();
    }

    public function measure_all(): static
    {
        if ($this->numClbits < $this->numQubits) {
            $this->numClbits = $this->numQubits;
        }
        // 既存の測定を残しつつ未測定のqubitを追加する
        for ($q = 0; $q < $this->numQubits; $q++) {
            if (!array_key_exists($q, $this->measureMap)) {
                $this->measureMap[$q] = $q;
                $this->operations[] = [
                    'name' => 'measure',
                    'qubits' => [$q],
                    'params' => [],
                    'clbits' => [$q],
                ];
            }
        }
        return $this;
    }

    public function barrier(int|array|\Traversable ...$qubits): static
    {
        $flat = [];
        foreach ($qubits as $q) {
            if ($q instanceof \Traversable) {
                $q = iterator_to_array($q, false);
            }
            if (is_array($q)) {
                foreach (array_values($q) as $v) {
                    $flat[] = (int)$v;
                }
            } else {
                $flat[] = (int)$q;
            }
        }
        if (count($flat) === 0) {
            $flat = range(0, $this->numQubits - 1);
        }
        foreach ($flat as $q) {
            $this->checkQubit($q);
        }
        // シミュレーション上は no-op として記録のみ
        $this->operations[] = [
            'name' => 'barrier',
            'qubits' => array_values($flat),
            'params' => [],
            'clbits' => [],
        ];
        return $this;
    }

    // ---- 情報・シミュレーション ----

    /** @return array<string,int> ゲート名 => 出現回数 (barrier除く) */
    public function countOps(): array
    {
        $counts = [];
        foreach ($this->operations as $op) {
            if ($op['name'] === 'barrier') {
                continue;
            }
            $counts[$op['name']] = ($counts[$op['name']] ?? 0) + 1;
        }
        return $counts;
    }

    public function depth(): int
    {
        // 簡易depth: qubitごとのレイヤー数を数える
        $layers = array_fill(0, $this->numQubits, 0);
        $maxDepth = 0;
        foreach ($this->operations as $op) {
            if ($op['name'] === 'barrier' || $op['name'] === 'measure') {
                continue;
            }
            $d = 0;
            foreach ($op['qubits'] as $q) {
                $d = max($d, $layers[$q]);
            }
            $d++;
            foreach ($op['qubits'] as $q) {
                $layers[$q] = $d;
            }
            $maxDepth = max($maxDepth, $d);
        }
        return $maxDepth;
    }

    /**
     * 最終状態ベクトルを返す. 各要素は ['re'=>float,'im'=>float].
     * @return array<int,array{re:float,im:float}>
     */
    public function statevector(): array
    {
        $sim = new StatevectorSimulator($this->numQubits);
        $sim->applyCircuit($this);
        return $sim->getStatevector();
    }

    /**
     * 計算基底の確率分布を返す. キーは qiskit形式ビット列 ('00' 等).
     * @return array<string,float>
     */
    public function probabilities(): array
    {
        $sim = new StatevectorSimulator($this->numQubits);
        $sim->applyCircuit($this);
        return $sim->probabilities($this->effectiveMeasureMap());
    }

    /**
     * @return array<int,int> qubit => clbit
     */
    public function effectiveMeasureMap(): array
    {
        if (!empty($this->measureMap)) {
            return $this->measureMap;
        }
        // 測定命令が無い場合は全qubit測定とみなす
        $map = [];
        for ($q = 0; $q < $this->numQubits; $q++) {
            $map[$q] = $q;
        }
        return $map;
    }

    public function effectiveNumClbits(): int
    {
        if ($this->numClbits > 0) {
            return $this->numClbits;
        }
        return $this->numQubits;
    }

    /** qiskitの draw() 的な簡易テキスト描画. 独立な操作は同一列に詰めて全ビットで幅を揃える */
    public function draw(): string
    {
        // 操作を列(レイヤー)にパッキング: 直前の列とqubitが重ならなければ同列に配置する.
        // これにより h([0,1]) や measure([0,1],[0,1]) が縦に揃って描画される.
        /** @var array<int,array{cells:array<int,string>,used:array<int,bool>,usedClbits:array<int,bool>,barrier:bool}> $columns */
        $columns = [];
        foreach ($this->operations as $op) {
            $label = $op['name'];
            if (in_array($op['name'], ['rx', 'ry', 'rz', 'p'], true)) {
                $label .= sprintf('(%.3g)', $op['params'][0]);
            } elseif ($op['name'] === 'u') {
                $label .= sprintf('(%.3g,%.3g,%.3g)', $op['params'][0], $op['params'][1], $op['params'][2]);
            }
            $qs = $op['qubits'];
            // この操作が使うセルと占有範囲を求める
            $contrib = []; // qubit => symbol
            $span = [];    // 占有するqubit範囲(2qubit間の│含む)
            $clbits = [];
            if ($op['name'] === 'measure') {
                $contrib[$qs[0]] = "M(c{$op['clbits'][0]})";
                $span = [$qs[0]];
                $clbits = [$op['clbits'][0]];
            } elseif ($op['name'] === 'barrier') {
                foreach ($qs as $q) {
                    $contrib[$q] = '‖';
                }
                $span = array_values($qs);
            } elseif (count($qs) === 1) {
                $contrib[$qs[0]] = $label;
                $span = [$qs[0]];
            } elseif ($op['name'] === 'cx' || $op['name'] === 'cnot') {
                [$c, $t] = $qs;
                $contrib[$c] = '●';
                $contrib[$t] = '⊕';
                $span = range(min($c, $t), max($c, $t));
                for ($q = min($c, $t) + 1; $q < max($c, $t); $q++) {
                    $contrib[$q] = '│';
                }
            } elseif ($op['name'] === 'cz') {
                foreach ($qs as $q) {
                    $contrib[$q] = '●';
                }
                [$a, $b] = $qs;
                $span = range(min($a, $b), max($a, $b));
                for ($q = min($a, $b) + 1; $q < max($a, $b); $q++) {
                    $contrib[$q] = '│';
                }
            } elseif ($op['name'] === 'swap') {
                foreach ($qs as $q) {
                    $contrib[$q] = '×';
                }
                [$a, $b] = $qs;
                $span = range(min($a, $b), max($a, $b));
                for ($q = min($a, $b) + 1; $q < max($a, $b); $q++) {
                    $contrib[$q] = '│';
                }
            } else {
                foreach ($qs as $q) {
                    $contrib[$q] = $label;
                }
                $span = array_values($qs);
            }

            $isBarrier = ($op['name'] === 'barrier');
            $merged = false;
            if (!$isBarrier && !empty($columns)) {
                $lastIdx = count($columns) - 1;
                if (!$columns[$lastIdx]['barrier']) {
                    $conflict = false;
                    foreach ($span as $q) {
                        if ($columns[$lastIdx]['used'][$q]) {
                            $conflict = true;
                            break;
                        }
                    }
                    if (!$conflict) {
                        foreach ($clbits as $c) {
                            if ($columns[$lastIdx]['usedClbits'][$c] ?? false) {
                                $conflict = true;
                                break;
                            }
                        }
                    }
                    if (!$conflict) {
                        foreach ($contrib as $q => $sym) {
                            $columns[$lastIdx]['cells'][$q] = $sym;
                        }
                        foreach ($span as $q) {
                            $columns[$lastIdx]['used'][$q] = true;
                        }
                        foreach ($clbits as $c) {
                            $columns[$lastIdx]['usedClbits'][$c] = true;
                        }
                        $merged = true;
                    }
                }
            }
            if (!$merged) {
                $cells = array_fill(0, $this->numQubits, '');
                foreach ($contrib as $q => $sym) {
                    $cells[$q] = $sym;
                }
                $used = array_fill(0, $this->numQubits, false);
                foreach ($span as $q) {
                    $used[$q] = true;
                }
                $usedClbits = [];
                foreach ($clbits as $c) {
                    $usedClbits[$c] = true;
                }
                $columns[] = ['cells' => $cells, 'used' => $used, 'usedClbits' => $usedClbits, 'barrier' => $isBarrier];
            }
        }

        $lines = [];
        for ($q = 0; $q < $this->numQubits; $q++) {
            $lines[$q] = "q{$q}: ─";
        }
        foreach ($columns as $col) {
            $cells = $col['cells'];
            // 列幅はこの列の最大表示幅. アイドルはその幅ぶんの '─' で埋める
            $width = 1;
            foreach ($cells as $cell) {
                if ($cell !== '') {
                    $width = max($width, self::displayWidth($cell));
                }
            }
            for ($q = 0; $q < $this->numQubits; $q++) {
                $cell = $cells[$q];
                if ($cell === '') {
                    $lines[$q] .= str_repeat('─', $width) . '─';
                } else {
                    $lines[$q] .= $cell . str_repeat('─', $width - self::displayWidth($cell)) . '─';
                }
            }
        }
        return implode("\n", $lines) . "\n";
    }

    private static function displayWidth(string $s): int
    {
        if (function_exists('mb_strwidth')) {
            return max(1, mb_strwidth($s, 'UTF-8'));
        }
        if (function_exists('mb_strlen')) {
            return max(1, mb_strlen($s, 'UTF-8'));
        }
        $chars = preg_split('//u', $s, -1, PREG_SPLIT_NO_EMPTY);
        if ($chars !== false) {
            return max(1, count($chars));
        }
        return max(1, strlen($s));
    }

    public function __toString(): string
    {
        return $this->draw();
    }
}
