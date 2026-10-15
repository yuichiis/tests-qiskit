<?php
declare(strict_types=1);

namespace Rindow\Qubit\Quantum;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;
use Rindow\Math\Matrix\MatrixOperator;
use Rindow\Qubit\QuantumCircuit;
use function Rindow\Math\Matrix\C;

/**
 * Rindow NDArray + BLAS 対応の状態ベクトルシミュレータ.
 * qubit 0 = LSB (Qiskit規約) で状態を保持する.
 *
 * 内部表現は complex128 NDArray x1: shape [dim]。
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
    private MatrixOperator $mo;
    private object $la;
    private NDArray $state;
    /** @var array<string,NDArray> 固定ゲート行列のデバイス側キャッシュ(同一$mo/同一デバイス内で再利用) */
    private array $gateCache = [];

    /** @var ?MatrixOperator サービス選択コスト回避のための共有インスタンス */
    private static ?MatrixOperator $sharedMo = null;

    public static function setSharedMatrixOperator(?MatrixOperator $mo): void
    {
        self::$sharedMo = $mo;
    }

    public static function sharedMatrixOperator(): MatrixOperator
    {
        if (self::$sharedMo === null) {
            self::$sharedMo = new MatrixOperator();
        }
        return self::$sharedMo;
    }

    public function __construct(int $numQubits, ?MatrixOperator $mo = null)
    {
        if ($numQubits < 1 || $numQubits > 30) {
            throw new InvalidArgumentException('numQubits must be in [1, 30].');
        }
        $this->numQubits = $numQubits;
        $this->dim = 1 << $numQubits;
        $this->mo = $mo ?? self::sharedMatrixOperator();
        $this->la = $this->mo->la();
        $this->state = $this->mo->zeros([$this->dim], NDArray::complex128);
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

    public function matrixOperator(): MatrixOperator
    {
        return $this->mo;
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

    /**
     * @param array<int> $qubits
     * @param array<float>|NDArray $params 回路由来はfloat配列。NDArray渡しはunitaryのみ対応(デバイス常駐用)。
     */
    public function applyOperation(string $name, array $qubits, array|NDArray $params = []): void
    {
        switch (strtolower($name)) {
            case 'h':
                $s = 1.0 / sqrt(2.0);
                $this->applySingle($qubits[0], $this->cachedGate('h', [[$s, $s], [$s, -$s]]));
                break;
            case 'x':
                $this->applySingle($qubits[0], $this->cachedGate('x', [[0, 1], [1, 0]]));
                break;
            case 'y':
                // [[0,-i],[i,0]]
                $this->applySingle($qubits[0], $this->cachedGate('y', [[0, C(0, -1)], [C(0, 1), 0]]));
                break;
            case 'z':
                $this->applySingle($qubits[0], $this->cachedGate('z', [[1, 0], [0, -1]]));
                break;
            case 's':
                $this->applySingle($qubits[0], $this->cachedGate('s', [[1, 0], [0, C(0, 1)]]));
                break;
            case 'sdg':
            case 's_dg':
                $this->applySingle($qubits[0], $this->cachedGate('sdg', [[1, 0], [0, C(0, -1)]]));
                break;
            case 't':
                $s = 1.0 / sqrt(2.0);
                $this->applySingle($qubits[0], $this->cachedGate('t', [[1, 0], [0, C($s, $s)]]));
                break;
            case 'tdg':
            case 't_dg':
                $s = 1.0 / sqrt(2.0);
                $this->applySingle($qubits[0], $this->cachedGate('tdg', [[1, 0], [0, C($s, -$s)]]));
                break;
            case 'id':
            case 'i':
                break;
            case 'rx':
                $this->applyRx($qubits[0], (float)$params[0]);
                break;
            case 'ry':
                $this->applyRy($qubits[0], (float)$params[0]);
                break;
            case 'rz':
                $this->applyRz($qubits[0], (float)$params[0]);
                break;
            case 'p':
                $this->applyPhase((float)$params[0], $qubits[0]);
                break;
            case 'u':
                $this->applyU((float)$params[0], (float)$params[1], (float)$params[2], $qubits[0]);
                break;
            case 'unitary':
                // NDArray直渡し(デバイス常駐ゲート用): shapeで1qubit/Nqubitを振り分け
                if ($params instanceof NDArray) {
                    $shape = $params->shape();
                    if ($shape === [2, 2]) {
                        $this->applySingle($qubits[0], $params);
                    } else {
                        $this->applyNQubit(array_reverse($qubits), $params);
                    }
                    break;
                }
                if (count($params) === 8) {
                    $this->applySingle($qubits[0], [
                        [C((float)$params[0], (float)$params[1]), C((float)$params[2], (float)$params[3])],
                        [C((float)$params[4], (float)$params[5]), C((float)$params[6], (float)$params[7])],
                    ]);
                    break;
                }
                // Nqubitユニタリ: paramsを行優先で復元し、qiskit順(LSB first)を
                // MSB firstに直して applyNQubit に渡す.
                // 回路側はfloatのまま、ここで一発だけNDArray化(H2Dは1回)する。
                $dim = (int)sqrt(count($params) / 2);
                $mat = [];
                $idx = 0;
                for ($i = 0; $i < $dim; $i++) {
                    $row = [];
                    for ($j = 0; $j < $dim; $j++) {
                        $row[] = C((float)$params[$idx], (float)$params[$idx + 1]);
                        $idx += 2;
                    }
                    $mat[] = $row;
                }
                $this->applyNQubit(array_reverse($qubits), $mat);
                break;
            case 'cx':
            case 'cnot':
                $this->applyCx($qubits[0], $qubits[1]);
                break;
            case 'cz':
                $this->applyCz($qubits[0], $qubits[1]);
                break;
            case 'cp':
                $this->applyCp((float)$params[0], $qubits[0], $qubits[1]);
                break;
            case 'swap':
                $this->applySwap($qubits[0], $qubits[1]);
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
     * 固定ゲート行列のデバイス側キャッシュ. 同一シミュレータ(同一$mo/デバイス)内では
     * $mo->array()によるH2D転送を初回1回だけにする.
     * @param array<array<float|object>> $matrix
     */
    private function cachedGate(string $key, array $matrix): NDArray
    {
        if (!isset($this->gateCache[$key])) {
            $this->gateCache[$key] = $this->mo->array($matrix, dtype: NDArray::complex128);
        }
        return $this->gateCache[$key];
    }

    /**
     * PHP行列 or デバイス常駐NDArrayをcomplex128 NDArrayに正規化する.
     * NDArrayの場合はH2D転送なしでそのまま使う(OpenCL対応の要点).
     * @param array<array<float|object>>|NDArray $matrix
     * @param array<int> $shape 期待shape ([2,2]/[4,4]/[dim,dim])
     */
    private function toUnitaryNDArray(array|NDArray $matrix, array $shape): NDArray
    {
        if ($matrix instanceof NDArray) {
            if ($matrix->shape() !== $shape) {
                throw new InvalidArgumentException(
                    'Matrix shape mismatch: [' . implode(',', $matrix->shape()) . '] vs [' . implode(',', $shape) . '].'
                );
            }
            if ($matrix->dtype() !== NDArray::complex128) {
                return $this->mo->astype($matrix, NDArray::complex128);
            }
            return $matrix;
        }
        return $this->mo->array($matrix, dtype: NDArray::complex128);
    }

    /**
     * 1qubitゲート U=[[a,b],[c,d]] を適用 (各要素は実数またはComplex、またはcomplex128 NDArray[2,2]).
     * @param array<array<float|object>>|NDArray $matrix 2x2 (row-major)
     */
    public function applySingle(int $qubit, array|NDArray $matrix): void
    {
        if ($qubit < 0 || $qubit >= $this->numQubits) {
            throw new InvalidArgumentException("Qubit index out of range: {$qubit}");
        }
        if ($matrix instanceof NDArray) {
            if ($matrix->shape() !== [2, 2]) {
                throw new InvalidArgumentException('applySingle() requires a 2x2 matrix.');
            }
        } elseif (count($matrix) !== 2 || count((array)$matrix[0]) !== 2 || count((array)$matrix[1]) !== 2) {
            throw new InvalidArgumentException('applySingle() requires a 2x2 matrix.');
        }
        $la = $this->la;
        $L = 1 << $qubit;
        $H = $this->dim >> ($qubit + 1);
        $M = $H * $L;
        $U = $this->toUnitaryNDArray($matrix, [2, 2]);

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
     * @param array<array<float|object>>|NDArray $mat 4x4 (row-major。各要素は実数またはComplex、またはcomplex128 NDArray[4,4])
     */
    private function applyTwoQubit(int $qb, int $qa, array|NDArray $mat): void
    {
        [$L, $M, $H] = $this->twoQubitDims($qb, $qa);
        $la = $this->la;
        $U = $this->toUnitaryNDArray($mat, [4, 4]);
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
            $this->applyTwoQubit($control, $target, $this->cachedGate('cx_hi',
                [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]));
        } else {
            // control=低位(qa): qb反転はqa=1内。00->00,01->11,10->10,11->01
            $this->applyTwoQubit($target, $control, $this->cachedGate('cx_lo',
                [[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]]));
        }
    }

    public function applyCz(int $q1, int $q2): void
    {
        $this->applyTwoQubit(max($q1, $q2), min($q1, $q2), $this->cachedGate('cz',
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]]));
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
        $this->applyTwoQubit(max($q1, $q2), min($q1, $q2), $this->cachedGate('swap',
            [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]));
    }

    /**
     * 任意のk qubitゲート (2^k x 2^k 複素行列) を適用する.
     * 基底順は $qubitsMsbFirst の先頭が最上位ビット (例: [qb,qa] で idx=qb*2+qa)。
     * 状態を行ベクトルと見て Y = X @ U^T で単一gemm実行する。
     * @param array<int> $qubitsMsbFirst 先頭がMSB (重複なし)
     * @param array<array<float|object>>|NDArray $matrix 2^k x 2^k (row-major。各要素は実数またはComplex、またはcomplex128 NDArray)
     */
    public function applyNQubit(array $qubitsMsbFirst, array|NDArray $matrix): void
    {
        $k = count($qubitsMsbFirst);
        if ($k < 1) {
            throw new InvalidArgumentException('applyNQubit() requires at least 1 qubit.');
        }
        $dim = 1 << $k;
        if ($matrix instanceof NDArray) {
            if ($matrix->shape() !== [$dim, $dim]) {
                throw new InvalidArgumentException('applyNQubit() matrix size mismatch.');
            }
        } else {
            if (count($matrix) !== $dim) {
                throw new InvalidArgumentException('applyNQubit() matrix size mismatch.');
            }
            foreach ($matrix as $row) {
                if (count((array)$row) !== $dim) {
                    throw new InvalidArgumentException('applyNQubit() matrix size mismatch.');
                }
            }
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
        $U = $this->toUnitaryNDArray($matrix, [$dim, $dim]);
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

    /**
     * デバイス常駐の状態ベクトルを返す(OpenCL時はデバイスメモリのまま)。
     * 数値データをホストのPHP array/floatに落とさない経路の正規出力。
     */
    public function getStatevectorNDArray(): NDArray
    {
        return $this->state;
    }

    /**
     * デバイス常駐の計算基底確率を返す。abs/squareはデバイス上で実行され、
     * D2H転送は呼び出し側がtoArray()するまで発生しない。
     */
    public function getProbabilitiesNDArray(): NDArray
    {
        // test.php の流儀: probs = square(abs(psi))
        return $this->la->square($this->la->abs($this->state));
    }

    /** @return array<int,array{re:float,im:float}> 後方互換: D2Hはこの1回(toArray)だけ */
    public function getStatevector(): array
    {
        $arr = $this->state->toArray();
        $out = [];
        foreach ($arr as $c) {
            if (is_array($c)) {
                // 万が一 [re,im] 形式で返るバックエンド向けのフォールバック
                $a = array_values($c);
                $out[] = ['re' => (float)($a[0] ?? 0.0), 'im' => (float)($a[1] ?? 0.0)];
            } else {
                $out[] = ['re' => (float)$c->real, 'im' => (float)$c->imag];
            }
        }
        return $out;
    }

    /** @return array<float> 計算基底の確率。D2HはtoArray()の1回だけ */
    public function getProbabilities(): array
    {
        return $this->getProbabilitiesNDArray()->toArray();
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
        $counts = [];
        for ($s = 0; $s < $this->dim; $s++) {
            $p = $stateProbs[$s];
            if ($p == 0.0) {
                continue;
            }
            $bits = array_fill(0, $numClbits, 0);
            foreach ($measureMap as $q => $c) {
                $bits[$c] = ($s >> $q) & 1;
            }
            $key = '';
            for ($c = $numClbits - 1; $c >= 0; $c--) {
                $key .= (string)$bits[$c];
            }
            $counts[$key] = ($counts[$key] ?? 0.0) + $p;
        }
        ksort($counts);
        return $counts;
    }

    /**
     * 確率分布からショットをサンプリングする.
     * Rindowの randomCategorical (Cレベル: cumsum+searchsorted+randomUniform) を使用し、
     * 大規模dim/shotsでも高速に動作する。test.php の流儀に沿う。
     * @param array<int,int> $measureMap
     * @return array<string,int> ビット列 => カウント
     */
    public function sampleCounts(array $measureMap, int $numClbits, int $shots, ?int $seed = null): array
    {
        if ($shots < 1) {
            throw new InvalidArgumentException('shots must be >= 1.');
        }
        // test.php の流儀: probs = square(abs(psi))
        $probsND = $this->la->scal(2.0, $this->la->log($this->la->abs($this->state)));
        $samplesND = $this->la->randomCategorical($probsND, numSamples: $shots, seed: $seed);
        $samplesArr = $samplesND->toArray();
        $counts = [];
        foreach ($samplesArr as $idx) {
            $idx = (int)$idx;
            $bits = array_fill(0, $numClbits, 0);
            foreach ($measureMap as $q => $c) {
                $bits[$c] = ($idx >> $q) & 1;
            }
            $key = '';
            for ($c = $numClbits - 1; $c >= 0; $c--) {
                $key .= (string)$bits[$c];
            }
            $counts[$key] = ($counts[$key] ?? 0) + 1;
        }
        ksort($counts);
        return $counts;
    }
}
