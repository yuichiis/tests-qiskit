<?php
declare(strict_types=1);

namespace RindowTest\Qubit;

use Rindow\Qubit\CircuitGate;

trait TestUtils
{
    /**
     * Assert a circuit gate matrix (row-major) against expected values.
     * Expected elements may be numbers, objects with real/imag, or
     * arrays with real/imag keys.
     *
     * @param array<array<mixed>> $expected
     */
    protected function assertGateMatrix(CircuitGate $gate, array $expected, float $delta = 1.0e-6): void
    {
        $expectedFlat = [];
        foreach ($expected as $row) {
            self::assertIsArray($row);
            foreach ($row as $value) {
                $expectedFlat[] = self::complexParts($value);
            }
        }

        $actual = $gate->parameters();
        self::assertCount(count($expectedFlat), $actual);
        foreach ($expectedFlat as $i => [$real, $imag]) {
            self::assertEqualsWithDelta(
                self::roundComplex($real),
                self::roundComplex($actual[$i]->real),
                $delta,
                "Real part mismatch at flat index {$i}."
            );
            self::assertEqualsWithDelta(
                self::roundComplex($imag),
                self::roundComplex($actual[$i]->imag),
                $delta,
                "Imaginary part mismatch at flat index {$i}."
            );
        }
    }

    /**
     * @param array<int,float> $expected
     * @param array<int,float> $actual
     */
    protected function assertFloatVectorEquals(array $expected, array $actual, float $delta = 1.0e-6): void
    {
        self::assertCount(count($expected), $actual);
        foreach ($expected as $i => $value) {
            self::assertEqualsWithDelta($value, $actual[$i], $delta, "Mismatch at index {$i}.");
        }
    }

    /**
     * @param array<int,array{float,float}> $expected [real, imag] pairs
     * @param array<int,object{real:float,imag:float}> $actual
     */
    protected function assertAmplitudesEqual(array $expected, array $actual, float $delta = 1.0e-6): void
    {
        self::assertCount(count($expected), $actual);
        foreach ($expected as $i => [$real, $imag]) {
            self::assertEqualsWithDelta($real, $actual[$i]->real, $delta, "Real amplitude mismatch at index {$i}.");
            self::assertEqualsWithDelta($imag, $actual[$i]->imag, $delta, "Imag amplitude mismatch at index {$i}.");
        }
    }

    /** @return array{float,float} */
    private static function complexParts(mixed $value): array
    {
        if (is_int($value) || is_float($value)) {
            return [(float)$value, 0.0];
        }
        if (is_object($value)) {
            if (isset($value->real)) {
                return [
                    (float)$value->real,
                    (float)($value->imag ?? 0.0),
                ];
            }
            if (method_exists($value, 'getReal')) {
                return [(float)$value->getReal(), (float)$value->getImag()];
            }
        }
        if (is_array($value)) {
            return [
                (float)($value['real'] ?? 0.0),
                (float)($value['imag'] ?? 0.0),
            ];
        }
        self::fail('Unsupported expected complex value.');
    }

    private static function roundComplex(float $value): float
    {
        return abs($value) < 1.0e-12 ? 0.0 : $value;
    }
}
