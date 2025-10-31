# Contributing to scBridge

Thank you for your interest in contributing to scBridge!

## Development Setup

### Python Development

```bash
cd python
pip install -e ".[dev]"
```

### R Development

```bash
cd R
# Install in development mode
R CMD INSTALL .
```

## Making Changes

1. Fork the repository
2. Create a new branch for your feature or bugfix
3. Make your changes
4. Test your changes
5. Submit a pull request

## Code Style

- **Python**: Follow PEP 8
- **R**: Follow tidyverse style guide

## Testing

- Python: Run tests with `pytest`
- R: Run tests with `devtools::test()`

## Pull Request Process

1. Update documentation if needed
2. Ensure all tests pass
3. Update [CHANGELOG.md](CHANGELOG.md) with your changes
4. Submit PR with clear description of changes

## Questions?

Open an issue or discussion on GitHub.
