<?php
declare(strict_types=1);

namespace Rindow\Qubit\Math;

use InvalidArgumentException;
use Interop\Polite\Math\Matrix\NDArray;

/**
 * 量子計算で使う浮動小数点精度のグローバル設定.
 *
 * 既定は float64 + complex128。float64 非対応の OpenCL 環境などで
 * float32 を使いたい場合は、計算を始める前に一度だけ設定する。
 *
 * ```php
 * use Interop\Polite\Math\Matrix\NDArray;
 * use Rindow\Qubit\Math\Precision;
 *
 * Precision::set(NDArray::float32);          // float32 + complex64
 * Precision::set(NDArray::complex128);       // float64 + complex128 (既定)
 * Precision::set(NDArray::float32, NDArray::complex64);
 * ```
 *
 * 実数 dtype と複素 dtype は必ず対応する組で使い、計算途中で
 * 切り替えないこと。ゲート行列・状態ベクトル・演算子行列などは
 * すべてこの設定から dtype を取得する。
 */
final class Precision
{
    private const REAL_TO_COMPLEX = [
        NDArray::float32 => NDArray::complex64,
        NDArray::float64 => NDArray::complex128,
    ];

    private const COMPLEX_TO_REAL = [
        NDArray::complex64 => NDArray::float32,
        NDArray::complex128 => NDArray::float64,
    ];

    private static int $realDtype = NDArray::float64;
    private static int $complexDtype = NDArray::complex128;

    /**
     * 既定の精度を設定する.
     *
     * 引数が1つの場合は float32/float64/complex64/complex128 のいずれかを
     * 受け取り、対応する組を設定する。2つの場合は実数 dtype と複素 dtype
     * の組を明示する (対応しない組は拒否する)。
     */
    public static function set(int $dtype, ?int $complexDtype = null): void
    {
        if ($complexDtype === null) {
            if (isset(self::REAL_TO_COMPLEX[$dtype])) {
                self::$realDtype = $dtype;
                self::$complexDtype = self::REAL_TO_COMPLEX[$dtype];
                return;
            }
            if (isset(self::COMPLEX_TO_REAL[$dtype])) {
                self::$realDtype = self::COMPLEX_TO_REAL[$dtype];
                self::$complexDtype = $dtype;
                return;
            }
            throw new InvalidArgumentException(
                'Precision dtype must be float32, float64, complex64, or complex128.'
            );
        }

        if (!isset(self::REAL_TO_COMPLEX[$dtype])) {
            throw new InvalidArgumentException('The real dtype must be float32 or float64.');
        }
        if (self::REAL_TO_COMPLEX[$dtype] !== $complexDtype) {
            throw new InvalidArgumentException(
                'The real and complex dtypes must form a matching pair: '
                . 'float32 + complex64 or float64 + complex128.'
            );
        }
        self::$realDtype = $dtype;
        self::$complexDtype = $complexDtype;
    }

    /** float32 + complex64 に切り替える. */
    public static function setFloat32(): void
    {
        self::set(NDArray::float32, NDArray::complex64);
    }

    /** float64 + complex128 に切り替える. */
    public static function setFloat64(): void
    {
        self::set(NDArray::float64, NDArray::complex128);
    }

    /** 既定値 (float64 + complex128) に戻す. */
    public static function reset(): void
    {
        self::setFloat64();
    }

    public static function realDtype(): int
    {
        return self::$realDtype;
    }

    public static function complexDtype(): int
    {
        return self::$complexDtype;
    }
}
