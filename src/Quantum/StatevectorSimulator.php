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

    /** @param array<int> $qubits @param array<float> $params */
    public function applyOperation(string $name, array $qubits, array $params = []): void
    {
        switch (strtolower($name)) {
            case 'h':
                $s = 1.0 / sqrt(2.0);
                $this->applySingle($qubits[0], [[$s, $s], [$s, -$s]]);
                break;
            case 'x':
                $this->applySingle($qubits[0], [[0, 1], [1, 0]]);
                break;
            case 'y':
                // [[0,-i],[i,0]]
                $this->applySingle($qubits[0], [[0, C(0, -1)], [C(0, 1), 0]]);
                break;
            case 'z':
                $this->applySingle($qubits[0], [[1, 0], [0, -1]]);
                break;
            case 's':
                $this->applySingle($qubits[0], [[1, 0], [0, C(0, 1)]]);
                break;
            case 'sdg':
            case 's_dg':
                $this->applySingle($qubits[0], [[1, 0], [0, C(0, -1)]]);
                break;
            case 't':
                $s = 1.0 / sqrt(2.0);
                $this->applySingle($qubits[0], [[1, 0], [0, C($s, $s)]]);
                break;
            case 'tdg':
            case 't_dg':
                $s = 1.0 / sqrt(2.0);
                $this->applySingle($qubits[0], [[1, 0], [0, C($s, -$s)]]);
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
                $this->applySingle($qubits[0], [
                    [C((float)$params[0], (float)$params[1]), C((float)$params[2], (float)$params[3])],
                    [C((float)$params[4], (float)$params[5]), C((float)$params[6], (float)$params[7])],
                ]);
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
     * 1qubitゲート U=[[a,b],[c,d]] を適用 (各要素は実数またはComplex).
     * @param array<array<float|object>> $matrix 2x2 (row-major)
     */
    public function applySingle(int $qubit, array $matrix): void
    {
        if ($qubit < 0 || $qubit >= $this->numQubits) {
            throw new InvalidArgumentException("Qubit index out of range: {$qubit}");
        }
        if (count($matrix) !== 2 || count((array)$matrix[0]) !== 2 || count((array)$matrix[1]) !== 2) {
            throw new InvalidArgumentException('applySingle() requires a 2x2 matrix.');
        }
        $la = $this->la;
        $mo = $this->mo;
        $L = 1 << $qubit;
        $H = $this->dim >> ($qubit + 1);
        $M = $H * $L;
        $U = $mo->array($matrix, dtype: NDArray::complex128);

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
        $U = $mo->array($mat, dtype: NDArray::complex128);
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

    /** @return array<int,array{re:float,im:float}> */
    public function getStatevector(): array
    {
        $buf = $this->state->buffer();
        $off = $this->state->offset();
        $out = [];
        for ($i = 0; $i < $this->dim; $i++) {
            $c = $buf[$off + $i];
            $out[] = ['re' => (float)$c->real, 'im' => (float)$c->imag];
        }
        return $out;
    }

    /** @return array<float> 計算基底の確率 */
    public function getProbabilities(): array
    {
        // test.php の流儀: probs = square(abs(psi))
        $probs = $this->la->square($this->la->abs($this->state));
        return $probs->toArray();
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
