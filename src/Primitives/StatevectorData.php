<?php
declare(strict_types=1);

namespace Rindow\Qubit\Primitives;

/**
 * qiskitの Statevector に相当. `getStatevector()` で最終状態ベクトルを返す.
 */
class StatevectorData
{
    /** @var array<int,object{real:float,imag:float}> */
    private array $statevector;

    /** @param array<int,object{real:float,imag:float}> $statevector */
    public function __construct(array $statevector)
    {
        $this->statevector = $statevector;
    }

    /**
     * @return array<int,object{real:float,imag:float}> 各要素は real/imag を持つ複素数オブジェクト
     */
    public function getStatevector(): array
    {
        return $this->statevector;
    }
}
