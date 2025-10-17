.PHONY: help install install-dev test test-cov lint format clean pre-commit-install

help:
	@echo "Available commands:"
	@echo "  make install           - Install package"
	@echo "  make install-dev       - Install package with dev dependencies"
	@echo "  make test              - Run tests"
	@echo "  make test-cov          - Run tests with coverage"
	@echo "  make lint              - Run linters (ruff, black check, isort check)"
	@echo "  make format            - Auto-format code with black and isort"
	@echo "  make pre-commit-install - Install pre-commit hooks"
	@echo "  make clean             - Remove build artifacts"

install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"

test:
	pytest tests/ -v

test-cov:
	pytest tests/ --cov=drhip --cov-report=html --cov-report=term

lint:
	@echo "Running ruff..."
	ruff check drhip/ tests/
	@echo "\nChecking formatting with black..."
	black --check drhip/ tests/
	@echo "\nChecking import sorting with isort..."
	isort --check-only drhip/ tests/

format:
	@echo "Formatting with black..."
	black drhip/ tests/
	@echo "Sorting imports with isort..."
	isort drhip/ tests/
	@echo "Fixing with ruff..."
	ruff check --fix drhip/ tests/

pre-commit-install:
	pre-commit install
	@echo "Pre-commit hooks installed! They will run automatically on git commit."

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
