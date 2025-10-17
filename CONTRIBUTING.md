# Contributing to DRHIP

Thank you for contributing to DRHIP (Data Reduction for HyPhy with Inference Processing)!

## Development Setup

### 1. Install Development Dependencies

```bash
# Install package with dev dependencies
pip install -e ".[dev]"

# Or use the Makefile
make install-dev
```

### 2. Install Pre-commit Hooks

Pre-commit hooks will automatically run linting and formatting checks before each commit:

```bash
# Install hooks
pre-commit install

# Or use the Makefile
make pre-commit-install
```

To manually run pre-commit on all files:
```bash
pre-commit run --all-files
```

## Running Tests

```bash
# Run all tests
make test

# Run tests with coverage report
make test-cov

# Or use pytest directly
pytest tests/ -v
pytest tests/ --cov=drhip --cov-report=html
```

## Code Quality

### Linting

```bash
# Check code quality (ruff, black, isort)
make lint

# Auto-format code
make format
```

### Manual Linting Commands

```bash
# Run ruff
ruff check drhip/ tests/

# Check formatting with black
black --check drhip/ tests/

# Check import sorting
isort --check-only drhip/ tests/
```

### Auto-fixing Issues

```bash
# Format code with black
black drhip/ tests/

# Sort imports with isort
isort drhip/ tests/

# Auto-fix ruff issues
ruff check --fix drhip/ tests/
```

## Pre-commit Hooks

The pre-commit hooks will automatically run:
- Trailing whitespace removal
- End-of-file fixer
- YAML/JSON/TOML validation
- Black formatting
- isort import sorting
- Ruff linting

If pre-commit fails, fix the issues and commit again. Most formatting issues can be auto-fixed with `make format`.

## Testing Guidelines

### Writing Tests

- Place tests in the `tests/` directory
- Name test files with `test_` prefix
- Use descriptive test function names starting with `test_`
- Include docstrings explaining what each test validates

### Backward Compatibility

When adding features that depend on specific HyPhy versions:
1. Make new fields optional in validation
2. Add appropriate warnings when features are unavailable
3. Create test data for both old and new versions
4. Add comprehensive tests to prevent regression

See `tests/test_cfel_backward_compatibility.py` for an example.

## Pull Request Process

1. Create a feature branch from `main`
2. Make your changes
3. Run tests: `make test`
4. Run linting: `make lint` (or let pre-commit handle it)
5. Commit your changes (pre-commit hooks will run automatically)
6. Push to your fork and create a pull request
7. GitHub Actions will automatically run tests on your PR

## Continuous Integration

GitHub Actions will automatically:
- Run tests on Python 3.8, 3.9, 3.10, 3.11
- Test on Ubuntu and macOS
- Generate coverage reports
- Run linting checks

All checks must pass before merging.

## Questions?

If you have questions or need help, please open an issue on GitHub.
