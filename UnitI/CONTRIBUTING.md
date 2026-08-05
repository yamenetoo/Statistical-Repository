# Contributing to unit_transform_dist

Contributions are welcome. This project favours small, well-tested pull
requests over large ones.

## Development setup

```bash
git clone https://github.com/yamenetoo/unit_transform_dist.git
cd unit_transform_dist
pip install -e ".[dev]"
```

## Running tests

```bash
pytest -v
```

## Style

- Keep functions type-hinted (the package ships `py.typed`).
- New estimators should raise `unit_transform_dist.FittingError` on
  optimisation failure, not a bare `RuntimeError`.
- Any function that samples or optimises stochastically should accept a
  `random_state` argument.
- Prefer adding a regression test alongside any bug fix, named
  `test_<bug>_<behaviour>` and referencing the original issue in a
  docstring, as in `tests/test_bugfixes.py`.

## Building and checking the distribution

```bash
pip install build twine
python -m build
twine check dist/*
```

## Reporting issues

Please include a minimal reproducible example (base distribution name,
data-generating code, and the exact call that fails).
