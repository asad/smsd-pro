# Contributing

Thanks for your interest in improving SMSD Pro!

## Ways to contribute
- Bug reports and minimal reproductions
- Performance benchmarks and difficult pairs
- New SMARTS predicate helpers ($name and $(...))
- Documentation improvements

## Development setup
1. Create an environment and install:
   ```bash
   python -m venv .venv && source .venv/bin/activate
   pip install -e .[rdkit]
   pip install -r requirements-dev.txt  # if present
   ```
2. Run tests and style checks:
   ```bash
   pytest -q
   ```

## Pull requests
- Keep changes focused; separate refactors from features
- Include tests or bench entries for new behaviour
- Update docs where helpful

By contributing, you agree that your contributions will be licensed under the Apache-2.0 license.
