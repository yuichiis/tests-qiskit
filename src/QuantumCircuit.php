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

    public function h(int $qubit): static
    {
        return $this->append('h', [$qubit]);
    }

    public function x(int $qubit): static
    {
        return $this->append('x', [$qubit]);
    }

    public function y(int $qubit): static
    {
        return $this->append('y', [$qubit]);
    }

    public function z(int $qubit): static
    {
        return $this->append('z', [$qubit]);
    }

    public function s(int $qubit): static
    {
        return $this->append('s', [$qubit]);
    }

    public function sdg(int $qubit): static
    {
        return $this->append('sdg', [$qubit]);
    }

    public function t(int $qubit): static
    {
        return $this->append('t', [$qubit]);
    }

    public function tdg(int $qubit): static
    {
        return $this->append('tdg', [$qubit]);
    }

    public function id(int $qubit): static
    {
        return $this->append('id', [$qubit]);
    }

    public function rx(float $theta, int $qubit): static
    {
        return $this->append('rx', [$qubit], [$theta]);
    }

    public function ry(float $theta, int $qubit): static
    {
        return $this->append('ry', [$qubit], [$theta]);
    }

    public function rz(float $theta, int $qubit): static
    {
        return $this->append('rz', [$qubit], [$theta]);
    }

    /** Phase gate P(λ) == U1(λ) */
    public function p(float $lambda, int $qubit): static
    {
        return $this->append('p', [$qubit], [$lambda]);
    }

    public function u(float $theta, float $phi, float $lambda, int $qubit): static
    {
        return $this->append('u', [$qubit], [$theta, $phi, $lambda]);
    }

    // ---- 2量子ビットゲート ----

    public function cx(int $control, int $target): static
    {
        if ($control === $target) {
            throw new InvalidArgumentException('Control and target must differ.');
        }
        return $this->append('cx', [$control, $target]);
    }

    /** alias of cx() */
    public function cnot(int $control, int $target): static
    {
        return $this->cx($control, $target);
    }

    public function cz(int $control, int $target): static
    {
        if ($control === $target) {
            throw new InvalidArgumentException('Control and target must differ.');
        }
        return $this->append('cz', [$control, $target]);
    }

    public function swap(int $qubit1, int $qubit2): static
    {
        if ($qubit1 === $qubit2) {
            throw new InvalidArgumentException('Swap qubits must differ.');
        }
        return $this->append('swap', [$qubit1, $qubit2]);
    }

    /**
     * 任意の1qubitユニタリ行列 [[a,b],[c,d]] を適用する.
     * 各要素は float(実数) または [re, im] / ['re'=>..,'im'=>..] / Complex オブジェクト.
     */
    public function unitary(array $matrix, int $qubit): static
    {
        $flat = self::matrixToComplexPairs($matrix);
        if (count($flat) !== 4) {
            throw new InvalidArgumentException('unitary() requires a 2x2 matrix.');
        }
        return $this->append('unitary', [$qubit], $flat);
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

    public function measure(int $qubit, int $clbit): static
    {
        $this->checkQubit($qubit);
        $this->checkClbit($clbit);
        $this->measureMap[$qubit] = $clbit;
        return $this->append('measure', [$qubit], [], [$clbit]);
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

    public function barrier(int ...$qubits): static
    {
        if (count($qubits) === 0) {
            $qubits = range(0, $this->numQubits - 1);
        }
        foreach ($qubits as $q) {
            $this->checkQubit($q);
        }
        // シミュレーション上は no-op として記録のみ
        $this->operations[] = [
            'name' => 'barrier',
            'qubits' => array_values($qubits),
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

    /** qiskitの draw() 的な簡易テキスト描画 */
    public function draw(): string
    {
        $lines = [];
        for ($q = 0; $q < $this->numQubits; $q++) {
            $lines[$q] = "q{$q}: ─";
        }
        foreach ($this->operations as $op) {
            $label = $op['name'];
            if (in_array($op['name'], ['rx', 'ry', 'rz', 'p'], true)) {
                $label .= sprintf('(%.3g)', $op['params'][0]);
            } elseif ($op['name'] === 'u') {
                $label .= sprintf('(%.3g,%.3g,%.3g)', $op['params'][0], $op['params'][1], $op['params'][2]);
            }
            $qs = $op['qubits'];
            if ($op['name'] === 'measure') {
                $lines[$qs[0]] .= "M(c{$op['clbits'][0]})─";
                continue;
            }
            if ($op['name'] === 'barrier') {
                foreach ($qs as $q) {
                    $lines[$q] .= "‖─";
                }
                continue;
            }
            if (count($qs) === 1) {
                $lines[$qs[0]] .= "{$label}─";
            } elseif ($op['name'] === 'cx' || $op['name'] === 'cnot') {
                [$c, $t] = $qs;
                $lines[$c] .= "●─";
                $lines[$t] .= "⊕─";
            } elseif ($op['name'] === 'cz') {
                foreach ($qs as $q) {
                    $lines[$q] .= "●─";
                }
            } elseif ($op['name'] === 'swap') {
                foreach ($qs as $q) {
                    $lines[$q] .= "×─";
                }
            } else {
                foreach ($qs as $q) {
                    $lines[$q] .= "{$label}─";
                }
            }
        }
        return implode("\n", $lines) . "\n";
    }

    public function __toString(): string
    {
        return $this->draw();
    }
}
