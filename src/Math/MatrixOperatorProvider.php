<?php
declare(strict_types=1);

namespace Rindow\Qubit\Math;

use Rindow\Math\Matrix\MatrixOperator;

final class MatrixOperatorProvider
{
    private static ?MatrixOperator $shared = null;

    public static function set(?MatrixOperator $mo): void
    {
        self::$shared = $mo;
    }

    public static function get(): MatrixOperator
    {
        if (self::$shared === null) {
            self::$shared = new MatrixOperator();
        }
        return self::$shared;
    }
}
