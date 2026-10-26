## PHP
```bash
phpvenv .phpvenv
. .\Scripts\Activate.ps1
composer require rindow/rindow-math-matrix-matlibffi

php samples\test_qubit.php
```

## Python
```bash
uv venv qiskit-venv
cd qiskit-venv
. .\Scripts\Activate.ps1
uv pip install qiskit matplotlib
deactivate

uv run --active python samples/test.py
```

### GPU
```
$env:RINDOW_QUBIT_LINEAR_ALGEBRA = "clblast::GPU"
php samples\test_qubit.php
```

