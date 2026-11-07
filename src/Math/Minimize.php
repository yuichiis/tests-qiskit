<?php
declare(strict_types=1);

namespace Rindow\Qubit\Math;

use InvalidArgumentException;

/**
 * scipy.optimize.minimize の簡易版.
 *
 * ```php
 * $res = Minimize::minimize($cost, $theta0, method: 'COBYLA', options: [
 *     'maxiter' => 500,
 *     'tol' => 1e-8,
 * ]);
 * echo $res->fun;
 * ```
 */
final class Minimize
{
    /**
     * @param callable(array<int,float>):float $fun
     * @param array<int,float> $x0
     * @param array<string,mixed>|null $options
     * @param array<int,callable(array<int,float>):float> $constraints
     */
    public static function minimize(
        callable $fun,
        array $x0,
        string $method = 'COBYLA',
        ?array $options = null,
        array $constraints = [],
    ): OptimizeResult {
        $options ??= [];
        $method = strtoupper($method);

        if ($method === 'COBYLA') {
            $maxiter = (int)($options['maxiter'] ?? 1000);
            $tol = (float)($options['tol'] ?? 1.0e-8);
            $rhobeg = (float)($options['rhobeg'] ?? 1.0);
            [$x, $f, $nfev] = Cobyla::minimize(
                $fun,
                $x0,
                $constraints,
                maxiter: $maxiter,
                tol: $tol,
                rhobeg: $rhobeg,
            );
            return new OptimizeResult(
                $x,
                $f,
                $nfev,
                true,
                0,
                'Optimization terminated successfully.',
            );
        }

        if ($method === 'NELDER-MEAD' || $method === 'NELDER_MEAD') {
            $maxiter = (int)($options['maxiter'] ?? 200);
            $tol = (float)($options['tol'] ?? 1.0e-6);
            $step = (float)($options['step'] ?? 0.3);
            [$x, $f] = NelderMead::minimize(
                $fun,
                $x0,
                step: $step,
                maxIter: $maxiter,
                tol: $tol,
            );
            return new OptimizeResult(
                $x,
                $f,
                $maxiter,
                true,
                0,
                'Optimization terminated successfully.',
            );
        }

        throw new InvalidArgumentException("Unsupported optimization method: {$method}");
    }
}
