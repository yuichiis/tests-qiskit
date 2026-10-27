<?php
declare(strict_types=1);

namespace RindowTest\Qubit\Circuit;

use InvalidArgumentException;
use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\TestCase;
use Rindow\Qubit\Circuit\AbstractMatrixGate;
use Rindow\Qubit\Circuit\MCPhaseGate;
use Rindow\Qubit\Circuit\MCXGate;
use RindowTest\Qubit\TestUtils;

class MultiControlledGateTest extends TestCase
{
    use TestUtils;

    public function testMcxOneControlIsCnot(): void
    {
        $gate = new MCXGate(1);

        $this->assertInstanceOf(AbstractMatrixGate::class, $gate);
        $this->assertSame('mcx', $gate->name());
        $this->assertSame(1, $gate->numCtrlQubits());
        $this->assertSame(2, $gate->numQubits());
        $this->assertSame(1, $gate->ctrlState());
        $this->assertSame(4, $gate->dimension());
        $this->assertGateMatrix($gate, [
            [1, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0],
            [0, 1, 0, 0],
        ]);
    }

    public function testMcxTwoControls(): void
    {
        $gate = new MCXGate(2);

        $this->assertSame(2, $gate->numCtrlQubits());
        $this->assertSame(3, $gate->numQubits());
        $this->assertSame(3, $gate->ctrlState());
        $this->assertSame(8, $gate->dimension());
        $this->assertGateMatrix($gate, [
            [1, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1],
            [0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 1, 0, 0, 0, 0],
        ]);
    }

    public function testMcxCustomControlState(): void
    {
        $gate = new MCXGate(2, 1);

        $this->assertSame(1, $gate->ctrlState());
        $expected = $this->identity(8);
        // control bits are 01, so inputs 1 and 5 are swapped through the target bit.
        $expected[1][5] = 1;
        $expected[5][1] = 1;
        $expected[1][1] = 0;
        $expected[5][5] = 0;
        $this->assertGateMatrix($gate, $expected);
    }

    public function testMcxDraw(): void
    {
        $this->assertSame(
            ['cells' => [0 => '●', 1 => '●', 2 => '⊕'], 'span' => [0, 1, 2]],
            (new MCXGate(2))->draw([0, 1, 2])
        );
    }

    public function testMcPhaseOneControlMatchesCPhase(): void
    {
        $lambda = M_PI / 3.0;
        $gate = new MCPhaseGate($lambda, 1);

        $this->assertInstanceOf(AbstractMatrixGate::class, $gate);
        $this->assertSame('mcphase', $gate->name());
        $this->assertSame(1, $gate->numCtrlQubits());
        $this->assertSame(2, $gate->numQubits());
        $this->assertSame(1, $gate->ctrlState());
        $this->assertSame($lambda, $gate->lambda());
        $this->assertSame(4, $gate->dimension());
        $this->assertGateMatrix($gate, [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, ['real' => cos($lambda), 'imag' => sin($lambda)]],
        ]);
    }

    public function testMcPhaseTwoControls(): void
    {
        $lambda = 0.7;
        $gate = new MCPhaseGate($lambda, 2);

        $this->assertSame(2, $gate->numCtrlQubits());
        $this->assertSame(3, $gate->numQubits());
        $this->assertSame(3, $gate->ctrlState());
        $this->assertSame(8, $gate->dimension());

        $expected = $this->identity(8);
        $expected[7][7] = ['real' => cos($lambda), 'imag' => sin($lambda)];
        $this->assertGateMatrix($gate, $expected);
    }

    public function testMcPhaseCustomControlState(): void
    {
        $lambda = 0.7;
        $gate = new MCPhaseGate($lambda, 2, 1);

        $this->assertSame(1, $gate->ctrlState());
        $expected = $this->identity(8);
        $expected[5][5] = ['real' => cos($lambda), 'imag' => sin($lambda)];
        $this->assertGateMatrix($gate, $expected);
    }

    public function testControlStateAsBinaryString(): void
    {
        $this->assertSame(2, (new MCXGate(2, '10'))->ctrlState());
        $this->assertSame(2, (new MCPhaseGate(0.5, 2, '10'))->ctrlState());
        $this->assertSame(0, (new MCXGate(2, '00'))->ctrlState());
    }

    public function testMcPhaseDraw(): void
    {
        $this->assertSame(
            ['cells' => [0 => '●', 1 => '●', 2 => 'P'], 'span' => [0, 1, 2]],
            (new MCPhaseGate(0.5, 2))->draw([0, 1, 2])
        );
    }

    public function testMcPhaseInverse(): void
    {
        $gate = new MCPhaseGate(0.5, 2, 1);
        $inverse = $gate->inverse();

        $this->assertInstanceOf(MCPhaseGate::class, $inverse);
        $this->assertSame(-0.5, $inverse->lambda());
        $this->assertSame(2, $inverse->numCtrlQubits());
        $this->assertSame(1, $inverse->ctrlState());
    }

    #[DataProvider('invalidConstructorProvider')]
    public function testInvalidConstructorArguments(callable $factory): void
    {
        $this->expectException(InvalidArgumentException::class);
        $factory();
    }

    /** @return array<string,array{callable}> */
    public static function invalidConstructorProvider(): array
    {
        return [
            'mcx zero controls' => [static fn() => new MCXGate(0)],
            'mcx too many controls' => [static fn() => new MCXGate(30)],
            'mcx negative control state' => [static fn() => new MCXGate(1, -1)],
            'mcx control state too large' => [static fn() => new MCXGate(1, 2)],
            'mcx bad binary string length' => [static fn() => new MCXGate(1, '00')],
            'mcx bad binary string character' => [static fn() => new MCXGate(2, '1x')],
            'mcphase zero controls' => [static fn() => new MCPhaseGate(0.5, 0)],
            'mcphase too many controls' => [static fn() => new MCPhaseGate(0.5, 30)],
            'mcphase negative control state' => [static fn() => new MCPhaseGate(0.5, 1, -1)],
            'mcphase control state too large' => [static fn() => new MCPhaseGate(0.5, 1, 2)],
            'mcphase bad binary string length' => [static fn() => new MCPhaseGate(0.5, 1, '00')],
            'mcphase bad binary string character' => [static fn() => new MCPhaseGate(0.5, 2, '1x')],
        ];
    }

    /** @return array<array<int>> */
    private function identity(int $dim): array
    {
        $matrix = [];
        for ($i = 0; $i < $dim; $i++) {
            $matrix[$i] = array_fill(0, $dim, 0);
            $matrix[$i][$i] = 1;
        }
        return $matrix;
    }
}
