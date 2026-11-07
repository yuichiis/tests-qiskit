<?php
declare(strict_types=1);

namespace Rindow\Qubit\Math;

/**
 * scipy.optimize.OptimizeResult の簡易版.
 */
class OptimizeResult
{
    /**
     * @param array<int,float> $x
     */
    public function __construct(
        public array $x,
        public float $fun,
        public int $nfev,
        public bool $success,
        public int $status,
        public string $message,
    ) {
    }
}
