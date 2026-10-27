<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Visualization;

use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Visualization\Visualization;

class VisualizationTest extends TestCase
{
    public function testEmptyCounts(): void
    {
        $this->assertSame("(empty counts)\n", Visualization::plotHistogram([]));
    }

    public function testHistogramIsSortedAndScaled(): void
    {
        $output = Visualization::plotHistogram(['1' => 1, '0' => 3], 10);
        $lines = explode("\n", rtrim($output, "\n"));

        $this->assertCount(2, $lines);
        $this->assertStringStartsWith('0 : ', $lines[0]);
        $this->assertStringStartsWith('1 : ', $lines[1]);
        $this->assertStringContainsString('    3', $lines[0]);
        $this->assertStringContainsString(' 75.0%', $lines[0]);
        $this->assertStringContainsString(str_repeat('█', 10), $lines[0]);
        $this->assertStringContainsString('    1', $lines[1]);
        $this->assertStringContainsString(' 25.0%', $lines[1]);
        $this->assertStringContainsString(str_repeat('█', 3), $lines[1]);
        $this->assertStringNotContainsString(str_repeat('█', 4), $lines[1]);
    }

    public function testAllZeroCounts(): void
    {
        $output = Visualization::plotHistogram(['0' => 0, '1' => 0], 10);
        $lines = explode("\n", rtrim($output, "\n"));

        $this->assertCount(2, $lines);
        $this->assertStringContainsString('  0.0%', $lines[0]);
        $this->assertStringNotContainsString('█', $output);
    }

    public function testShowHistogram(): void
    {
        $this->expectOutputString(Visualization::plotHistogram(['0' => 1], 4));
        Visualization::showHistogram(['0' => 1], 4);
    }
}
