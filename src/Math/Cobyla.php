<?php
declare(strict_types=1);

namespace Rindow\Qubit\Math;

use InvalidArgumentException;

/**
 * COBYLA (Constrained Optimization BY Linear Approximation) 風の
 * 導関数不要オプティマイザ.
 *
 * 現在の最良点の周りに単体 (simplex) を作り、その線形補間から
 * 勾配を推定して信頼領域内で降下する。制約を渡した場合はペナルティ
 * 法で扱うため、厳密な COBYLA の線形計画部分は簡略化している。
 * VQE のような滑らかな無制約問題を想定している。
 *
 * ```php
 * $result = Cobyla::minimize($cost, $theta0, maxiter: 500, tol: 1e-8);
 * ```
 */
final class Cobyla
{
    /**
     * @param callable(array<int,float>):float $fun
     * @param array<int,float> $x0
     * @param array<int,callable(array<int,float>):float> $constraints 各制約は >= 0 で満たす
     * @return array{array<int,float>,float,int} [最良点, 目的関数値, 関数評価回数]
     */
    public static function minimize(
        callable $fun,
        array $x0,
        array $constraints = [],
        int $maxiter = 1000,
        float $tol = 1.0e-8,
        float $rhobeg = 1.0,
    ): array {
        if (count($x0) === 0) {
            throw new InvalidArgumentException('x0 must not be empty.');
        }
        if ($maxiter < 0) {
            throw new InvalidArgumentException('maxiter must be >= 0.');
        }
        if ($tol <= 0.0) {
            throw new InvalidArgumentException('tol must be > 0.');
        }
        if ($rhobeg <= 0.0) {
            throw new InvalidArgumentException('rhobeg must be > 0.');
        }

        $n = count($x0);
        $x = [];
        foreach ($x0 as $value) {
            $x[] = (float)$value;
        }

        $nfev = 0;
        $evaluate = static function (array $point) use ($fun, $constraints, &$nfev): array {
            $nfev++;
            $raw = (float)$fun($point);
            $merit = $raw;
            foreach ($constraints as $constraint) {
                $violation = max(0.0, -(float)$constraint($point));
                $merit += 1.0e6 * $violation * $violation;
            }
            return [$raw, $merit];
        };

        [$fx, $merit] = $evaluate($x);
        $bestX = $x;
        $bestFun = $fx;
        $bestMerit = $merit;

        $rho = $rhobeg;
        $rhoEnd = max($tol, 1.0e-14);
        // 勾配の差分幅。rho が小さくなりすぎて差分が数値誤差に埋もれるのを防ぐ。
        $minGradientStep = 1.0e-4;

        for ($iter = 0; $iter < $maxiter && $rho >= $rhoEnd; $iter++) {
            $h = max($rho, $minGradientStep);
            $gradient = [];
            $norm2 = 0.0;
            for ($i = 0; $i < $n; $i++) {
                $probe = $x;
                $probe[$i] += $h;
                [, $probeMerit] = $evaluate($probe);
                $g = ($probeMerit - $merit) / $h;
                $gradient[] = $g;
                $norm2 += $g * $g;
            }
            $norm = sqrt($norm2);
            if ($norm < 1.0e-14) {
                // 一次のモデルが平坦。これ以上改善できる方向が無い。
                $rho *= 0.5;
                continue;
            }

            $candidate = [];
            foreach ($gradient as $i => $g) {
                $candidate[$i] = $x[$i] - ($rho / $norm) * $g;
            }
            [$candidateFun, $candidateMerit] = $evaluate($candidate);

            if ($candidateMerit < $merit) {
                $x = $candidate;
                $fx = $candidateFun;
                $merit = $candidateMerit;
                if ($merit < $bestMerit) {
                    $bestMerit = $merit;
                    $bestFun = $fx;
                    $bestX = $x;
                }
                // 成功したら信頼領域を少し広げて加速する (上限は rhobeg)。
                $rho = min($rhobeg, $rho * 1.5);
            } else {
                $rho *= 0.5;
            }
        }

        return [$bestX, $bestFun, $nfev];
    }
}
