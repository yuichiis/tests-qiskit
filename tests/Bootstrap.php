<?php
declare(strict_types=1);

date_default_timezone_set('UTC');

$projectVendor = __DIR__ . '/../vendor/autoload.php';
if (file_exists($projectVendor)) {
    $loader = require $projectVendor;
} else {
    if (!getenv('COMPOSER_HOME')) {
        // When PHPUnit is installed globally, derive the Composer home from
        // PHPUnit's own location so the test suite can run without extra env vars.
        if (class_exists(\PHPUnit\Framework\TestCase::class)) {
            $phpunitFile = (new \ReflectionClass(\PHPUnit\Framework\TestCase::class))->getFileName();
            if ($phpunitFile !== false) {
                $candidate = dirname($phpunitFile, 6);
                if (file_exists($candidate . '/vendor/autoload.php')) {
                    putenv('COMPOSER_HOME=' . $candidate);
                }
            }
        }
    }
    $loader = require __DIR__ . '/init_autoloader.php';
}

$loader->addPsr4('RindowTest\\Qubit\\', __DIR__ . '/RindowTest/Qubit');

return $loader;
