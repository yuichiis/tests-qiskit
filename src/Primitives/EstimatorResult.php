<?php
declare(strict_types=1);

namespace Rindow\Qubit\Primitives;

use ArrayAccess;
use Countable;
use IteratorAggregate;
use LogicException;
use OutOfBoundsException;
use Rindow\Qubit\QuantumInfo\Statevector;
use Traversable;

/**
 * qiskit の PrimitiveResult に相当. `$result[0]->data->evs` で期待値を取得する.
 */
class EstimatorResult implements ArrayAccess, Countable, IteratorAggregate
{
    /** @var array<int,EstimatorPubResult> */
    private array $pubs;

    /** @param array<int,EstimatorPubResult> $pubs */
    public function __construct(array $pubs)
    {
        $this->pubs = array_values($pubs);
    }

    public function get(int $index): EstimatorPubResult
    {
        if (!isset($this->pubs[$index])) {
            throw new OutOfBoundsException("Pub index out of range: {$index}");
        }
        return $this->pubs[$index];
    }

    public function getStatevector(int $index = 0): Statevector
    {
        return $this->get($index)->getStatevector();
    }

    public function getExpectationValue(int $index = 0): float
    {
        return $this->get($index)->getExpectationValue();
    }

    // --- ArrayAccess / Countable / IteratorAggregate ---

    public function offsetExists(mixed $offset): bool
    {
        return isset($this->pubs[$offset]);
    }

    public function offsetGet(mixed $offset): EstimatorPubResult
    {
        return $this->get((int)$offset);
    }

    public function offsetSet(mixed $offset, mixed $value): void
    {
        throw new LogicException('EstimatorResult is immutable.');
    }

    public function offsetUnset(mixed $offset): void
    {
        throw new LogicException('EstimatorResult is immutable.');
    }

    public function count(): int
    {
        return count($this->pubs);
    }

    public function getIterator(): Traversable
    {
        return new \ArrayIterator($this->pubs);
    }
}
