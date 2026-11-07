<?php
declare(strict_types=1);

namespace Rindow\Qubit\Circuit;

use InvalidArgumentException;
use Rindow\Qubit\QuantumCircuit;

/**
 * qiskit.circuit.library.real_amplitudes のPHP版.
 *
 * RY ゲートだけを使う実数振幅 ansatz。パラメータ数は (reps + 1) * numQubits。
 *
 * ```php
 * $ansatz = new RealAmplitudes(numQubits: 2, reps: 2, entanglement: 'linear');
 * $theta = array_fill(0, $ansatz->numParameters(), 0.0);
 * $circuit = $ansatz->assignParameters($theta);
 * ```
 */
class RealAmplitudes
{
    public const ENTANGLEMENT_FULL = 'full';
    public const ENTANGLEMENT_LINEAR = 'linear';
    public const ENTANGLEMENT_REVERSE_LINEAR = 'reverse_linear';
    public const ENTANGLEMENT_PAIRWISE = 'pairwise';
    public const ENTANGLEMENT_CIRCULAR = 'circular';

    private int $numQubits;
    private int $reps;
    /** @var array<int,array{0:int,1:int}> */
    private array $entanglementPairs;
    private string $parameterPrefix;
    private ?object $la;

    /**
     * @param string|array<int,array{0:int,1:int}> $entanglement
     */
    public function __construct(
        int $numQubits,
        int $reps = 3,
        string|array $entanglement = self::ENTANGLEMENT_FULL,
        string $parameterPrefix = 'θ',
        ?object $la = null
    ) {
        if ($numQubits < 1) {
            throw new InvalidArgumentException('numQubits must be >= 1.');
        }
        if ($reps < 0) {
            throw new InvalidArgumentException('reps must be >= 0.');
        }
        $this->numQubits = $numQubits;
        $this->reps = $reps;
        $this->parameterPrefix = $parameterPrefix;
        $this->la = $la;
        $this->entanglementPairs = $this->buildEntanglement($numQubits, $entanglement);
    }

    /**
     * @param string|array<int,array{0:int,1:int}> $entanglement
     * @return array<int,array{0:int,1:int}>
     */
    private function buildEntanglement(int $numQubits, string|array $entanglement): array
    {
        if (is_array($entanglement)) {
            $pairs = [];
            foreach ($entanglement as $pair) {
                if (!is_array($pair) || count($pair) !== 2) {
                    throw new InvalidArgumentException('Entanglement pairs must be [control, target].');
                }
                [$a, $b] = array_values($pair);
                $a = (int)$a;
                $b = (int)$b;
                if ($a < 0 || $a >= $numQubits || $b < 0 || $b >= $numQubits || $a === $b) {
                    throw new InvalidArgumentException("Invalid entanglement pair: ({$a}, {$b})");
                }
                $pairs[] = [$a, $b];
            }
            return $pairs;
        }

        return match ($entanglement) {
            self::ENTANGLEMENT_FULL => $this->fullPairs($numQubits),
            self::ENTANGLEMENT_LINEAR => $this->linearPairs($numQubits),
            self::ENTANGLEMENT_REVERSE_LINEAR => array_map(
                static fn(array $pair): array => [$pair[1], $pair[0]],
                $this->linearPairs($numQubits)
            ),
            self::ENTANGLEMENT_PAIRWISE,
            self::ENTANGLEMENT_CIRCULAR => $this->pairwisePairs($numQubits),
            default => throw new InvalidArgumentException(
                "Unsupported entanglement: {$entanglement}"
            ),
        };
    }

    /** @return array<int,array{0:int,1:int}> */
    private function fullPairs(int $numQubits): array
    {
        $pairs = [];
        for ($i = 0; $i < $numQubits; $i++) {
            for ($j = $i + 1; $j < $numQubits; $j++) {
                $pairs[] = [$i, $j];
            }
        }
        return $pairs;
    }

    /** @return array<int,array{0:int,1:int}> */
    private function linearPairs(int $numQubits): array
    {
        $pairs = [];
        for ($i = 0; $i < $numQubits - 1; $i++) {
            $pairs[] = [$i, $i + 1];
        }
        return $pairs;
    }

    /** @return array<int,array{0:int,1:int}> */
    private function pairwisePairs(int $numQubits): array
    {
        $pairs = $this->linearPairs($numQubits);
        if ($numQubits > 2) {
            $pairs[] = [$numQubits - 1, 0];
        }
        return $pairs;
    }

    public function numQubits(): int
    {
        return $this->numQubits;
    }

    /** @deprecated use numQubits() */
    public function num_qubits(): int
    {
        return $this->numQubits;
    }

    public function reps(): int
    {
        return $this->reps;
    }

    /** @return array<int,array{0:int,1:int}> */
    public function entanglementPairs(): array
    {
        return $this->entanglementPairs;
    }

    public function numParameters(): int
    {
        return ($this->reps + 1) * $this->numQubits;
    }

    /** @deprecated use numParameters() */
    public function num_parameters(): int
    {
        return $this->numParameters();
    }

    /**
     * パラメータを記号 (θ[0] など) のまま保持した回路を返す。
     * 描画や構造確認用。シミュレーションには assignParameters() を使う。
     */
    public function decompose(): QuantumCircuit
    {
        $circuit = new QuantumCircuit($this->numQubits, la: $this->la);
        $parameter = 0;
        for ($rep = 0; $rep <= $this->reps; $rep++) {
            for ($q = 0; $q < $this->numQubits; $q++) {
                $label = $this->parameterPrefix . '[' . $parameter . ']';
                $circuit->append('ry', [$q], [$label]);
                $parameter++;
            }
            if ($rep < $this->reps) {
                foreach ($this->entanglementPairs as [$control, $target]) {
                    $circuit->cx($control, $target);
                }
            }
        }
        return $circuit;
    }

    /**
     * パラメータ値を割り当てた実行可能な回路を返す。
     *
     * @param array<int,float> $theta
     */
    public function assignParameters(array $theta): QuantumCircuit
    {
        $theta = array_values($theta);
        $expected = $this->numParameters();
        if (count($theta) !== $expected) {
            throw new InvalidArgumentException(
                'Parameter count mismatch: expected ' . $expected . ', got ' . count($theta) . '.'
            );
        }
        $circuit = new QuantumCircuit($this->numQubits, la: $this->la);
        $parameter = 0;
        for ($rep = 0; $rep <= $this->reps; $rep++) {
            for ($q = 0; $q < $this->numQubits; $q++) {
                $circuit->ry((float)$theta[$parameter], $q);
                $parameter++;
            }
            if ($rep < $this->reps) {
                foreach ($this->entanglementPairs as [$control, $target]) {
                    $circuit->cx($control, $target);
                }
            }
        }
        return $circuit;
    }

    /** assignParameters() の別名 */
    public function bind(array $theta): QuantumCircuit
    {
        return $this->assignParameters($theta);
    }

    public function draw(): string
    {
        return $this->decompose()->draw();
    }

    public function __toString(): string
    {
        return $this->draw();
    }
}
