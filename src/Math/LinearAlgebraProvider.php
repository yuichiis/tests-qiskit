<?php
declare(strict_types=1);

namespace Rindow\Qubit\Math;

use Rindow\Math\Matrix\MatrixOperator;
use Interop\Polite\Math\Matrix\OpenCL;

final class LinearAlgebraProvider
{
    private static ?object $shared = null;

    public static function set(?object $la): void
    {
        self::$shared = $la;
    }

    public static function get(): object
    {
        if (self::$shared === null) {
            $mo = new MatrixOperator();
            self::$shared = self::getLinearAlgebra($mo);
        }
        return self::$shared;
    }

    public static function getLinearAlgebra(
        MatrixOperator $mo,
        ?string $options = null
    ): object
    {
        if($options===null) {
            $env = getenv('RINDOW_QUBIT_LINEAR_ALGEBRA');
            if(empty($env)) {
                return $mo->la();
            } else {
                return self::getLinearAlgebra($mo,$env);
            }
        }
        
        $options = explode('::',$options);
        $serviceName = array_shift($options);
        if($serviceName!=='clblast') {
            throw new RuntimeException('Invalid service name');
        }
        $options = implode('::',$options);
        $options = strtoupper($options);
        if($options=='GPU') {
            $options = ['deviceType' => OpenCL::CL_DEVICE_TYPE_GPU];
        } elseif($options=='CPU') {
            $options = ['deviceType' => OpenCL::CL_DEVICE_TYPE_CPU];
        } elseif($options!=='') {
            $options = ['device' => $options];
        } else {
            $options = null;
        }
        return $mo->laAccelerated($serviceName,$options);
    }

}
