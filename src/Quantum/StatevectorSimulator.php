<?php
declare(strict_types=1);

namespace Rindow\Qubit\Quantum;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Math\Matrix\MatrixOperator;
use Rindow\Qubit\Circuit\CPhaseGate;
use Rindow\Qubit\Circuit\CXGate;
use Rindow\Qubit\Circuit\CZGate;
use Rindow\Qubit\Circuit\HGate;
use Rindow\Qubit\Circuit\IGate;
use Rindow\Qubit\Circuit\PhaseGate;
use Rindow\Qubit\Circuit\RXGate;
use Rindow\Qubit\Circuit\RYGate;
use Rindow\Qubit\Circuit\RZGate;
use Rindow\Qubit\Circuit\SdgGate;
use Rindow\Qubit\Circuit\SGate;
use Rindow\Qubit\Circuit\SwapGate;
use Rindow\Qubit\Circuit\TdgGate;
use Rindow\Qubit\Circuit\TGate;
use Rindow\Qubit\Circuit\UGate;
use Rindow\Qubit\Circuit\XGate;
use Rindow\Qubit\Circuit\YGate;
use Rindow\Qubit\Circuit\ZGate;
use Rindow\Qubit\CircuitGate;
use Rindow\Qubit\QuantumCircuit;
use Rindow\Qubit\Math\IntMath;
use Rindow\Qubit\Math\LinearAlgebraProvider;
use function Rindow\Math\Matrix\C;

/**
 * Rindow NDArray + BLAS 対応の状態ベクトルシミュレータ.
 * qubit 0 = LSB (Qiskit規約) で状態を保持する.
 *
 * 内部表現は complex64 NDArray x1: shape [dim]。
 * rindow-matlib の transpose が複素数の共役なし転置で
 * N次元対応になったため、1/2qubitゲートとも転置+gemm各1発で書ける。
 * test.php の流儀 (MatrixOperator + LinearAlgebra::gemm 等) に沿う。
 *
 * 1qubitゲート U(2x2):
 *   [H,2,L] -> transpose[0,2,1] -> [H,L,2] -> reshape[M,2] (M=H*L)
 *   Y = X @ UT (UT=転置のみ、共役なし、単一zgemm)
 *   -> reshape[H,L,2] -> transpose[0,2,1] -> [H,2,L] -> reshape[dim]
 * 2qubitゲート U(4x4、基底順 (qb,qa)):
 *   [H,2,M,2,L] -> transpose[0,2,4,1,3] -> [H,M,L,2,2] -> reshape[N,4]
 *   Y = X @ UT -> reshape -> transpose[0,3,1,4,2] -> reshape[dim]
 */
class StatevectorSimulator
{
    private int $numQubits;
    private int $dim;
    private object $la;
    private NDArray $state;

    public static function sharedLinearAlgebra(): object
    {
        return LinearAlgebraProvider::get();
    }

    public function __construct(
        int $numQubits,
        ?object $la = null
    ) {
        if ($numQubits < 1 || $numQubits > 30) {
            throw new InvalidArgumentException('numQubits must be in [1, 30].');
        }
        $this->numQubits = $numQubits;
        $this->dim = 1 << $numQubits;
        $this->la = $la ?? LinearAlgebraProvider::get();
        $this->state = $this->la->zeros([$this->dim], NDArray::complex64);
        $this->state[0] = C(1, 0);
    }

    public function numQubits(): int
    {
        return $this->numQubits;
    }

    public function dimension(): int
    {
        return $this->dim;
    }

    /** テスト・デバッグ用: 内部NDArrayへのアクセス */
    public function stateNDArray(): NDArray
    {
        return $this->state;
    }

    public function la(): object
    {
        return $this->la;
    }

    public function applyCircuit(QuantumCircuit $circuit): void
    {
        if ($circuit->numQubits() !== $this->numQubits) {
            throw new InvalidArgumentException('Circuit qubit count mismatch.');
        }
        foreach ($circuit->operations() as $op) {
            $this->applyOperation($op['name'], $op['qubits'], $op['params']);
        }
    }

    /** @param array<int> $qubits @param array<mixed>|NDArray $params */
    public function applyOperation(string $name, array $qubits, array|NDArray $params = []): void
    {
        $name = strtolower($name);
        if ($params instanceof NDArray) {
            if (in_array($name, ['measure', 'barrier', 'reset'], true)) {
                return;
            }
            $this->applyMatrix($qubits, $params);
            return;
        }
        $gate = $this->standardGate($name, $params);
        if ($gate !== null) {
            $this->applyMatrix($qubits, $gate->matrix());
            return;
        }
        switch ($name) {
            case 'unitary':
            case 'qft':
            case 'qft_dg':
            case 'mcx':
            case 'mcphase':
                // 直接呼び出し互換のため、旧形式の複素オブジェクト配列を受け付ける。
                $this->applyMatrix($qubits, $this->gateMatrixNDArray(
                    $this->unitaryParamsToMatrix($params)
                ));
                break;
            case 'measure':
            case 'barrier':
            case 'reset':
                break;
            default:
                throw new InvalidArgumentException("Unknown gate: {$name}");
        }
    }

    /**
     * @param array<int> $qubits Qiskit order (first argument = LSB).
     * @param array<mixed> $params
     */
    private function standardGate(string $name, array $params): ?CircuitGate
    {
        return match ($name) {
            'h' => new HGate($this->mo),
            'x' => new XGate($this->mo),
            'y' => new YGate($this->mo),
            'z' => new ZGate($this->mo),
            's' => new SGate($this->mo),
            'sdg', 's_dg' => new SdgGate($this->mo),
            't' => new TGate($this->mo),
            'tdg', 't_dg' => new TdgGate($this->mo),
            'id', 'i' => new IGate($this->mo),
            'rx' => new RXGate((float)$params[0], $this->mo),
            'ry' => new RYGate((float)$params[0], $this->mo),
            'rz' => new RZGate((float)$params[0], $this->mo),
            'p' => new PhaseGate((float)$params[0], $this->mo),
            'u' => new UGate((float)$params[0], (float)$params[1], (float)$params[2], $this->mo),
            'cx', 'cnot' => new CXGate($this->mo),
            'cz' => new CZGate($this->mo),
            'cp' => new CPhaseGate((float)$params[0], $this->mo),
            'swap' => new SwapGate($this->mo),
            default => null,
        };
    }

    /** @param array<int> $qubits Qiskit order (first argument = LSB). */
    private function applyMatrix(array $qubits, NDArray $matrix): void
    {
        $U = $this->gateMatrixNDArray($matrix);
        $dim = $U->shape()[0];
        $expected = 1 << count($qubits);
        if ($dim !== $expected) {
            throw new InvalidArgumentException(
                "Gate matrix is {$dim}x{$dim} but " . count($qubits) . ' qubit(s) were given.'
            );
        }
        if ($dim === 2) {
            $this->applySingle($qubits[0], $U);
            return;
        }
        $this->applyNQubit(array_reverse($qubits), $U);
    }

    /** @param array<array<float|object>> $params 旧形式 (行優先) */
    private function unitaryParamsToMatrix(array $params): array
    {
        $dim = (int)sqrt(count($params));
        if ($dim * $dim !== count($params)) {
            throw new InvalidArgumentException('unitary() params size mismatch.');
        }
        $mat = [];
        $idx = 0;
        for ($i = 0; $i < $dim; $i++) {
            $row = [];
            for ($j = 0; $j < $dim; $j++) {
                $c = $params[$idx];
                $row[] = C((float)$c->real, (float)$c->imag);
                $idx++;
            }
            $mat[] = $row;
        }
        return $mat;
    }

    private function gateMatrixNDArray(array|NDArray $matrix): NDArray
    {
        $U = $matrix instanceof NDArray
            ? $matrix
            : $this->la->array($matrix, dtype: NDArray::complex64);
        $shape = $U->shape();
        if (count($shape) !== 2 || $shape[0] !== $shape[1]) {
            throw new InvalidArgumentException('Gate matrix must be a square 2D matrix.');
        }
        if ($U->dtype() !== NDArray::complex64) {
            $U = $this->la->astype($U, NDArray::complex64);
        }
        return $U;
    }

    /**
     * 1qubitゲート U=[[a,b],[c,d]] を適用 (各要素は実数またはComplex).
     * @param array<array<float|object>>|NDArray $matrix 2x2 (row-major)
     */
    public function applySingle(int $qubit, array|NDArray $matrix): void
    {
        if ($qubit < 0 || $qubit >= $this->numQubits) {
            throw new InvalidArgumentException("Qubit index out of range: {$qubit}");
        }
        $U = $this->gateMatrixNDArray($matrix);
        if ($U->shape() !== [2, 2]) {
            throw new InvalidArgumentException('applySingle() requires a 2x2 matrix.');
        }
        $la = $this->la;
        $L = 1 << $qubit;
        $H = $this->dim >> ($qubit + 1);
        $M = $H * $L;

        // [H,2,L] view -> transpose -> [H,L,2] -> [M,2]
        $flat = $la->transpose($this->state->reshape([$H, 2, $L]), perm: [0, 2, 1])
            ->reshape([$M, 2]);
        // Y = X @ U^T (単一zgemm。transB=共役なし転置。conj省略時は共役転置になるので明示)
        $out = $la->gemm($flat, $U, transB: true, conjB: false);
        // [M,2] -> [H,L,2] -> transpose -> [H,2,L] -> [dim]
        $this->state = $la->transpose($out->reshape([$H, $L, 2]), perm: [0, 2, 1])
            ->reshape([$this->dim]);
    }

    public function applyRx(int $qubit, float $theta): void
    {
        $c = cos($theta / 2.0);
        $s = sin($theta / 2.0);
        // [[c, -i s], [-i s, c]]
        $this->applySingle($qubit, [[$c, C(0, -$s)], [C(0, -$s), $c]]);
    }

    public function applyRy(int $qubit, float $theta): void
    {
        $c = cos($theta / 2.0);
        $s = sin($theta / 2.0);
        $this->applySingle($qubit, [[$c, -$s], [$s, $c]]);
    }

    public function applyRz(int $qubit, float $theta): void
    {
        $c = cos($theta / 2.0);
        $s = sin($theta / 2.0);
        $this->applySingle($qubit, [[C($c, -$s), 0], [0, C($c, $s)]]);
    }

    public function applyPhase(float $lambda, int $qubit): void
    {
        $this->applySingle($qubit, [[1, 0], [0, C(cos($lambda), sin($lambda))]]);
    }

    public function applyU(float $theta, float $phi, float $lambda, int $qubit): void
    {
        $c = cos($theta / 2.0);
        $s = sin($theta / 2.0);
        // [[c, -e^{iλ}s],[e^{iφ}s, e^{i(φ+λ)}c]]
        $this->applySingle($qubit, [
            [$c, C(-$s * cos($lambda), -$s * sin($lambda))],
            [C($s * cos($phi), $s * sin($phi)), C($c * cos($phi + $lambda), $c * sin($phi + $lambda))],
        ]);
    }

    /**
     * 2qubitインデックス分解: a=min, b=max として
     * L=2^a, M=2^(b-a-1), H=2^(n-b-1)。
     * 状態 index = ((h*2+qb)*M + m)*2*L + qa*L + l。
     * @return array{int,int,int} [L,M,H]
     */
    private function twoQubitDims(int $q1, int $q2): array
    {
        if ($q1 < 0 || $q1 >= $this->numQubits || $q2 < 0 || $q2 >= $this->numQubits) {
            throw new InvalidArgumentException("Qubit index out of range: {$q1}, {$q2}");
        }
        if ($q1 === $q2) {
            throw new InvalidArgumentException('Qubits must differ.');
        }
        $a = min($q1, $q2);
        $b = max($q1, $q2);
        $L = 1 << $a;
        $M = 1 << ($b - $a - 1);
        $H = $this->dim >> ($b + 1);
        return [$L, $M, $H];
    }

    /**
     * 2qubitゲート (4x4複素行列) を適用。基底順は (qb,qa) = (高位,低位) で
     * idx = qb*2+qa (00,01,10,11)。行ベクトル形式 X @ U^T で単一zgemm実行。
     * @param array<array<float|object>> $mat 4x4 (row-major。各要素は実数またはComplex)
     */
    private function applyTwoQubit(int $qb, int $qa, array $mat): void
    {
        [$L, $M, $H] = $this->twoQubitDims($qb, $qa);
        $la = $this->la;
        $mo = $this->mo;
        $U = $mo->array($mat, dtype: NDArray::complex64);
        $N = $H * $M * $L;
        // [H,qb,M,qa,L] -> [H,M,L,qb,qa] -> [N,4]
        $flat = $la->transpose(
            $this->state->reshape([$H, 2, $M, 2, $L]),
            perm: [0, 2, 4, 1, 3]
        )->reshape([$N, 4]);
        // Y = X @ U^T (transB=共役なし転置。conj省略時は共役転置になるので明示)
        $out = $la->gemm($flat, $U, transB: true, conjB: false);
        // [N,4] -> [H,M,L,2,2] -> [H,2,M,2,L] -> [dim]
        $this->state = $la->transpose(
            $out->reshape([$H, $M, $L, 2, 2]),
            perm: [0, 3, 1, 4, 2]
        )->reshape([$this->dim]);
    }

    public function applyCx(int $control, int $target): void
    {
        if ($control === $target) {
            throw new InvalidArgumentException('Control and target must differ.');
        }
        if ($control > $target) {
            // control=高位(qb): qa反転はqb=1内。00->00,01->01,10->11,11->10
            $this->applyTwoQubit($control, $target,
                [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]);
        } else {
            // control=低位(qa): qb反転はqa=1内。00->00,01->11,10->10,11->01
            $this->applyTwoQubit($target, $control,
                [[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]]);
        }
    }

    public function applyCz(int $q1, int $q2): void
    {
        $this->applyTwoQubit(max($q1, $q2), min($q1, $q2),
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]]);
    }

    /**
     * 制御位相ゲート CP(θ) = diag(1,1,1,e^{iθ}).
     * 両方が1のときのみ位相が付くため control/target の順序に依らない.
     */
    public function applyCp(float $theta, int $q1, int $q2): void
    {
        if ($q1 === $q2) {
            throw new InvalidArgumentException('Control and target must differ.');
        }
        $this->applyTwoQubit(max($q1, $q2), min($q1, $q2),
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
             [0, 0, 0, C(cos($theta), sin($theta))]]);
    }

    public function applySwap(int $q1, int $q2): void
    {
        $this->applyTwoQubit(max($q1, $q2), min($q1, $q2),
            [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]);
    }

    /**
     * 任意のk qubitゲート (2^k x 2^k 複素行列) を適用する.
     * 基底順は $qubitsMsbFirst の先頭が最上位ビット (例: [qb,qa] で idx=qb*2+qa)。
     * 状態を行ベクトルと見て Y = X @ U^T で単一gemm実行する。
     * @param array<int> $qubitsMsbFirst 先頭がMSB (重複なし)
     * @param array<array<float|object>>|NDArray $matrix 2^k x 2^k (row-major。各要素は実数またはComplex)
     */
    public function applyNQubit(array $qubitsMsbFirst, array|NDArray $matrix): void
    {
        $k = count($qubitsMsbFirst);
        if ($k < 1) {
            throw new InvalidArgumentException('applyNQubit() requires at least 1 qubit.');
        }
        $dim = 1 << $k;
        $U = $this->gateMatrixNDArray($matrix);
        if ($U->shape() !== [$dim, $dim]) {
            throw new InvalidArgumentException('applyNQubit() matrix size mismatch.');
        }
        $seen = [];
        foreach ($qubitsMsbFirst as $q) {
            if ($q < 0 || $q >= $this->numQubits || isset($seen[$q])) {
                throw new InvalidArgumentException("Qubit index out of range or duplicated: {$q}");
            }
            $seen[$q] = true;
        }
        $la = $this->la;
        $n = $this->numQubits;
        // 状態を [2]*n テンソルと見る. axis a (左から0始まり) <-> qubit n-1-a.
        $toAxis = static fn(int $q): int => $n - 1 - $q;
        $targetAxes = array_map($toAxis, array_values($qubitsMsbFirst));
        $isTarget = array_fill(0, $n, false);
        foreach ($targetAxes as $a) {
            $isTarget[$a] = true;
        }
        $otherAxes = [];
        for ($a = 0; $a < $n; $a++) {
            if (!$isTarget[$a]) {
                $otherAxes[] = $a;
            }
        }
        // 対象axisを末尾に集める (MSB firstの順序を保つ)
        $perm = array_merge($otherAxes, $targetAxes);
        $shape = array_fill(0, $n, 2);
        $rest = 1 << ($n - $k);
        $flat = $la->transpose($this->state->reshape($shape), perm: $perm)
            ->reshape([$rest, $dim]);
        // Y = X @ U^T (transB=共役なし転置。conj省略時は共役転置になるので明示)
        $out = $la->gemm($flat, $U, transB: true, conjB: false);
        // 逆置換で元に戻す
        $backShape = array_merge(
            array_fill(0, $n - $k, 2),
            array_fill(0, $k, 2)
        );
        $inv = array_fill(0, $n, 0);
        foreach ($perm as $new => $old) {
            $inv[$old] = $new;
        }
        $this->state = $la->transpose($out->reshape($backShape), perm: $inv)
            ->reshape([$this->dim]);
    }

    /** @return array<int,object{real:float,imag:float}> */
    public function getStatevector(): array
    {
        $flat = $this->state->reshape([$this->state->size()]);
        $out = [];
        foreach ($flat as $c) {
            $value = $this->la->scalar($c);
            $out[] = (object)['real' => (float)$value->real, 'imag' => (float)$value->imag];
        }
        return $out;
    }

    /** @return NDArray 計算基底の確率 (float64, shape [dim]) */
    public function getProbabilities(): NDArray
    {
        // test.php の流儀: probs = square(abs(psi))
        return $this->la->square($this->la->abs($this->state));
    }

    /** @param array<int,int> $measureMap */
    private function isIdentityMeasureMap(array $measureMap): bool
    {
        if (count($measureMap) !== $this->numQubits) {
            return false;
        }
        for ($q = 0; $q < $this->numQubits; $q++) {
            if (!array_key_exists($q, $measureMap) || (int)$measureMap[$q] !== $q) {
                return false;
            }
        }
        return true;
    }

    /**
     * 測定マップを考慮したビット列 => 確率.
     * @param array<int,int> $measureMap qubit => clbit
     * @param int $numClbits
     * @return array<string,float>
     */
    public function probabilities(array $measureMap = [], int $numClbits = 0): array
    {
        $stateProbs = $this->getProbabilities();
        if (empty($measureMap)) {
            for ($q = 0; $q < $this->numQubits; $q++) {
                $measureMap[$q] = $q;
            }
        }
        if ($numClbits <= 0) {
            $numClbits = max(array_values($measureMap)) + 1;
        }
        // 恒等マップ (q -> q) なら state index がそのまま計測ビット列になるため、
        // index変換表と scatterAdd を省き、reshape した確率をそのまま読む。
        if ($this->isIdentityMeasureMap($measureMap) && $numClbits >= $this->numQubits) {
            $counts = [];
            foreach ($stateProbs->reshape([$this->dim])->toArray() as $s => $p) {
                if ($p == 0.0) {
                    continue;
                }
                $counts[IntMath::int2bit($s, $numClbits)] = $p;
            }
            return $counts;
        }
        // state index -> clbit index の変換表を作り、scatterAdd で集計する.
        // measureMap に同じ clbit が複数回現れた場合は後勝ち (従来動作)。
        $indices = $this->la->zeros([$this->dim], NDArray::int32);
        for ($s = 0; $s < $this->dim; $s++) {
            $t = 0;
            foreach ($measureMap as $q => $c) {
                if ($c < 0 || $c >= $numClbits) {
                    continue;
                }
                $bit = ($s >> $q) & 1;
                $t = ($t & ~(1 << $c)) | ($bit << $c);
            }
            $indices[$s] = $t;
        }
        $scattered = $this->la->zeros([1 << $numClbits], NDArray::float64);
        $this->la->scatterAdd($indices, $stateProbs, $scattered);
        $counts = [];
        foreach ($scattered->toArray() as $t => $p) {
            if ($p == 0.0) {
                continue;
            }
            $counts[IntMath::int2bit($t, $numClbits)] = $p;
        }
        return $counts;
    }

    /**
     * 確率分布からショットをサンプリングする.
     * Rindowの randomCategorical (Cレベル: cumsum+searchsorted+randomUniform) でサンプリングし、
     * カウントは bincount (C/OpenCLレベル) で集計する. ホスト側では
     * 非ゼロ状態のビット列キー変換だけを行う (test.php の流儀)。
     * @param array<int,int> $measureMap
     * @return array<string,int> ビット列 => カウント
     */
    public function sampleCounts(array $measureMap, int $numClbits, int $shots, ?int $seed = null): array
    {
        if ($shots < 1) {
            throw new InvalidArgumentException('shots must be >= 1.');
        }
        // test.php の流儀: probs = square(abs(psi))
        $probsLogitsND = $this->la->scal(2.0, $this->la->log($this->la->abs($this->state)));
        $samplesND = $this->la->randomCategorical($probsLogitsND, numSamples: $shots, seed: $seed);
        // 状態indexごとの出現件数を bincount で集計 (shots個のPHPカウントループ不要)
        $stateCounts = $this->la->bincount($samplesND, $this->dim)->toArray();
        $counts = [];
        foreach ($stateCounts as $s => $n) {
            if ($n == 0) {
                continue;
            }
            $t = 0;
            foreach ($measureMap as $q => $c) {
                if ($c < 0 || $c >= $numClbits) {
                    continue;
                }
                $bit = ((int)$s >> $q) & 1;
                $t = ($t & ~(1 << $c)) | ($bit << $c);
            }
            $key = IntMath::int2bit($t, $numClbits);
            $counts[$key] = ($counts[$key] ?? 0) + (int)$n;
        }
        ksort($counts);
        return $counts;
    }
}
